function [theta_hat, r_hat, p_hat, u_hat, info] = nf_clkl(R_hat, W, P)
%NF_CLKL  Curvature-Learned KL Covariance Fitting estimator (proposed).
%
%  IMPROVEMENTS in v5 (from paper audits):
%  _________________________________________
%  1. MULTI-START WARM-START (3 runs, from Zhang et al. insight):
%     At high SNR, the OMP warm-start locks into a local minimum ~40% of
%     trials (N_MC=200 diagnostics: conv%=14% at SNR=10, 57% conv but
%     -5.03 dB at SNR=20 vs -5.65 dB at SNR=10 -- convergence to wrong min).
%     Fix: run OMP warm-start at 3 different u-initialisation points
%     (u_near, u_mid, u_far), evaluate the KL objective for each,
%     and keep the best support. O(3) overhead on warm-start only.
%
%  2. PER-ATOM RING-INDEXED u INITIALISATION (from Cui & Dai insight):
%     Instead of starting all active atoms at u_mid, each atom at angle
%     theta_k is initialised at the nearest Fresnel coherence ring:
%       u_k0 = s_k / (Z__ * sin^2(theta_k)),
%     where s_k = argmin_s | (1/s)*Z__*sin^2theta_k - r_mid |
%     and Z__ = D^2/(2beta^2lambda) is the Fresnel threshold distance.
%     For our small-array regime (Z__ < r_min), s_k = 1 for all theta,
%     so this reduces to u_k0 = 1/(Z__*sin^2theta_k) clipped to [u_min,u_max].
%     This gives curvature-aware initialisation relative to each atom's
%     angle, unlike the flat u_mid used in v4.
%
%  Changelog (cumulative):
%  v1: N0 frozen (eigenvalue estimate)
%  v2: tol_clkl relaxed; Q_theta grid increased
%  v3: OMP warm start with covariance deflation
%  v4: power-only main loop + 4-pass global residual MF scan
%  v5 (this): multi-start warm-start + per-atom ring-indexed init

% ---- Unpack ---------------------------------------------------------
M        = P.M;
N_RF     = P.N_RF;
Q_theta  = P.Q_theta;
d        = P.d;
lambda   = P.lambda;
d_ant    = P.d_ant;
lam_reg  = P.lambda_reg;
max_iter = P.max_iter;
tol      = P.tol_clkl;
alpha_p0 = P.alpha_p / N_RF;

% ---- Ablation flags (all default false; do NOT set in production runs) ---
% P.ablation_update_N0    -- re-estimate N0 every 10 iters (v1 collapse mode)
% P.ablation_single_start -- skip starts 1-2, use only Start 3 (far-range)
% P.ablation_skip_scan    -- skip Phase 2 post-loop MF scan
abl_N0     = isfield(P,'ablation_update_N0')    && P.ablation_update_N0;
abl_single = isfield(P,'ablation_single_start') && P.ablation_single_start;
abl_noscan = isfield(P,'ablation_skip_scan')    && P.ablation_skip_scan;

D_ap  = (M-1)*d_ant;
r_RD  = 2*D_ap^2/lambda;
r_min = P.r_lo_fac * r_RD / P.u_margin;
r_max = P.r_hi_fac * r_RD * P.u_margin;
u_min = 1/r_max;
u_max = 1/r_min;

theta_grid = linspace(P.theta_lo, P.theta_hi, Q_theta).';
c_vec      = (pi*d_ant^2/lambda) * sin(theta_grid).^2;
m_bar      = ((0:M-1).' - (M-1)/2);
m_bar2     = m_bar.^2;
omega_vec  = (2*pi*d_ant/lambda) * cos(theta_grid).';
lin_phase  = m_bar * omega_vec;    % M x Q_theta
I_NRF      = eye(N_RF);

R_hat_sym = (R_hat + R_hat') / 2;

% ====================================================================
%  N0: eigenvalue estimate, frozen
% ====================================================================
ev_sorted  = sort(real(eig(R_hat_sym)), 'ascend');
n_noise_ev = max(1, N_RF - d);
N0_hat     = max(mean(ev_sorted(1:n_noise_ev)), 1e-12);

% ====================================================================
%  Per-atom ring-indexed u initialisation (v5, Cui & Dai insight)
%  ---------------------------------------------------------------
%  Fresnel threshold distance Z__ = D^2/(2beta^2lambda).
%  For each grid angle theta_n, the ring-1 curvature init is:
%    u_ring(theta_n) = 1 / (Z__ * sin^2(theta_n))  clipped to [u_min, u_max]
%  When Z__*sin^2theta < r_min (our small-array regime), u_ring > u_max,
%  so the clip sets u_ring(theta_n) = u_max for near-boresight angles.
% ====================================================================
beta_delta = P.beta_delta;
Z_delta    = D_ap^2 / (2 * beta_delta^2 * lambda);

u_ring = zeros(Q_theta, 1);
for n = 1:Q_theta
    sin2_n = sin(theta_grid(n))^2;
    if sin2_n > 1e-6
        u_ring(n) = min(u_max, max(u_min, 1/(Z_delta * sin2_n)));
    else
        u_ring(n) = u_min;   % near-boresight: far-field init
    end
end

u_mid   = 0.5 * (u_min + u_max);
u_near  = u_max;    % near-range initialisation
u_far   = u_min;    % far-range initialisation

% ====================================================================
%  Multi-start OMP warm-start (v5)
%  Three u-initialisations: ring-indexed, near, far
%  Keep the one that achieves the lowest KL objective.
% ====================================================================
u_inits = {u_ring, u_near*ones(Q_theta,1), u_far*ones(Q_theta,1)};
best_active = [];
best_p_init = [];
best_L      = Inf;

% Ablation: single-start skips ring/near inits, uses only Start 3 (far-range)
si_range = 1:3;
if abl_single; si_range = 3; end

for si = si_range
    u_init_si = u_inits{si};
    D_si      = build_D(lin_phase, c_vec, u_init_si, m_bar2, W, N_RF, Q_theta);

    R_res_si  = R_hat_sym;
    active_si = zeros(d, 1);
    p_si      = zeros(Q_theta, 1);

    for k = 1:d
        sc_k = real(sum(conj(D_si) .* (R_res_si * D_si), 1));
        [~, best_k] = max(sc_k);
        active_si(k) = best_k;

        d_k   = D_si(:, best_k);
        nd2_k = max(real(d_k'*d_k), 1e-15);
        p_k   = max(0, (real(d_k'*R_res_si*d_k) - N0_hat*nd2_k)) / nd2_k^2;
        p_si(best_k) = p_k;

        Pk       = d_k * d_k' / nd2_k;
        R_res_si = R_res_si - Pk*R_res_si - R_res_si*Pk + Pk*R_res_si*Pk;
        R_res_si = (R_res_si + R_res_si') / 2;
    end

    % Evaluate KL objective for this warm-start
    Ry_si = D_si*diag(p_si)*D_si' + N0_hat*I_NRF;
    Ry_si = (Ry_si + Ry_si')/2;
    L_si  = kl_obj(Ry_si, R_hat_sym, p_si, lam_reg);

    if L_si < best_L
        best_L      = L_si;
        best_active = active_si;
        best_p_init = p_si;
        best_u_init = u_init_si;
    end
end

p_vec = best_p_init;
p_sc  = max(p_vec);
if p_sc > 0; p_vec = p_vec / p_sc; else; p_vec = zeros(Q_theta,1); end

% Build D_mat at the best warm-start's u values
D_mat = build_D(lin_phase, c_vec, best_u_init, m_bar2, W, N_RF, Q_theta);

% ====================================================================
%  Main loop: POWER UPDATES ONLY (D_mat fixed)
% ====================================================================
Ry     = D_mat*diag(p_vec)*D_mat' + N0_hat*I_NRF;
Ry     = (Ry+Ry')/2;
L_prev = kl_obj(Ry, R_hat_sym, p_vec, lam_reg);
L_hist = nan(max_iter, 1);
converged = false;

for t = 1:max_iter
    Ry_reg = Ry + 1e-10*I_NRF;
    Ry_inv = inv(Ry_reg);   %#ok<MINV>
    G      = Ry_inv - Ry_inv*R_hat_sym*Ry_inv;

    GD         = G * D_mat;
    grad_p_vec = real(sum(conj(D_mat).*GD, 1)).' + lam_reg;

    alpha_p = alpha_p0;
    descent = grad_p_vec.'*grad_p_vec;
    for ls = 1:12
        p_try  = max(0, p_vec - alpha_p*grad_p_vec);
        Ry_try = D_mat*diag(p_try)*D_mat' + N0_hat*I_NRF;
        Ry_try = (Ry_try+Ry_try')/2;
        L_try  = kl_obj(Ry_try, R_hat_sym, p_try, lam_reg);
        if L_try <= L_prev - P.ls_sigma*alpha_p*descent; break; end
        alpha_p = alpha_p*P.ls_beta;
    end
    p_vec = p_try;
    Ry    = Ry_try;

    % Ablation: re-estimate N0 from residual eigenvalues every 10 iters.
    % This is the documented v1 behaviour that causes gradient collapse at
    % high SNR (N0_hat drifts to floor ~1e-4 as p_vec dominates residual).
    if abl_N0 && mod(t, 10) == 0
        R_res_abl  = (R_hat_sym - D_mat*diag(p_vec)*D_mat');
        R_res_abl  = (R_res_abl + R_res_abl') / 2;
        ev_res_abl = sort(real(eig(R_res_abl)), 'ascend');
        N0_hat     = max(mean(ev_res_abl(1:n_noise_ev)), 1e-12);
        Ry         = D_mat*diag(p_vec)*D_mat' + N0_hat*I_NRF;
        Ry         = (Ry + Ry') / 2;
        L_prev     = kl_obj(Ry, R_hat_sym, p_vec, lam_reg); % sync cost
    end

    L_curr     = kl_obj(Ry, R_hat_sym, p_vec, lam_reg);
    L_hist(t)  = L_curr;
    rel_change = abs(L_curr - L_prev) / (abs(L_prev) + 1e-15);
    L_prev     = L_curr;
    if rel_change < tol && t > 5
        converged = true;
        break
    end
end

p_hat = p_vec;
[~, idx_sorted] = sort(p_hat, 'descend');
active_final    = idx_sorted(1:d);

% ====================================================================
%  POST-LOOP: 4-pass global 2D residual MF scan (v4, unchanged)
%  Skipped when abl_noscan=true; warm-start u values used directly.
% ====================================================================
Q_scan_th = 192;
Q_scan_u  = 256;

theta_ref = theta_grid(active_final);
u_ref     = best_u_init(active_final);   % v5: use ring-indexed init, not u_mid

d_ref = zeros(N_RF, d);
for k = 1:d
    d_ref(:,k) = atom(theta_ref(k), u_ref(k), m_bar, m_bar2, d_ant, lambda, W);
end

if ~abl_noscan
th_scan_global = linspace(P.theta_lo, P.theta_hi, Q_scan_th);
u_scan_global  = linspace(u_min, u_max, Q_scan_u);

for pass = 1:4
    for k = 1:d
        R_k = R_hat_sym;
        for j = 1:d
            if j == k; continue; end
            dj  = d_ref(:,j);
            R_k = R_k - p_hat(active_final(j)) * (dj*dj');
        end
        R_k = (R_k + R_k')/2;

        if mod(pass,2) == 1
            u_k = u_ref(k);
            c_k = (pi*d_ant^2/lambda)*sin(th_scan_global).^2;
            kap = c_k .* u_k;
            PHI = m_bar*(2*pi*d_ant/lambda*cos(th_scan_global)) - m_bar2*kap;
            Ds  = W' * exp(1j*PHI);
            sc  = real(sum(conj(Ds).*(R_k*Ds), 1));
            [~,bi] = max(sc);
            theta_ref(k) = th_scan_global(bi);
        else
            th_k    = theta_ref(k);
            omega_k = (2*pi*d_ant/lambda)*cos(th_k);
            c_k_sc  = (pi*d_ant^2/lambda)*sin(th_k)^2;
            kap     = c_k_sc * u_scan_global;
            PHI     = omega_k*m_bar - m_bar2*kap;
            Ds      = W' * exp(1j*PHI);
            sc      = real(sum(conj(Ds).*(R_k*Ds), 1));
            [~,bi]  = max(sc);
            u_ref(k) = u_scan_global(bi);
        end

        d_ref(:,k) = atom(theta_ref(k), u_ref(k), m_bar, m_bar2, d_ant, lambda, W);
    end
end
end  % ~abl_noscan

theta_hat = theta_ref;
r_hat     = 1 ./ u_ref;

u_hat_full = u_mid * ones(Q_theta, 1);
u_hat_full(active_final) = u_ref;
u_hat = u_hat_full;

info.L_hist      = L_hist(1:t);
info.n_iter      = t;
info.N0_hat      = N0_hat;
info.N0_init     = N0_hat;
info.active      = active_final;
info.converged   = converged;
info.p_init      = p_hat;
u_candidates = {u_ring, u_near*ones(Q_theta,1), u_far*ones(Q_theta,1)};
bs = 1;
for bsi = 1:3
    if isequal(u_candidates{bsi}, best_u_init); bs = bsi; break; end
end
info.best_start  = bs;
end


function d_k = atom(theta, u, m_bar, m_bar2, d_ant, lambda, W)
omega = (2*pi*d_ant/lambda)*cos(theta);
kappa = (pi*d_ant^2/lambda)*sin(theta)^2*u;
a     = exp(1j*omega*m_bar - 1j*kappa*m_bar2);
d_k   = W' * a;
end

function D = build_D(lin_phase, c_vec, u_vec, m_bar2, W, N_RF, Q)
kappa_vec = c_vec .* u_vec;
PHI       = lin_phase - m_bar2 * kappa_vec.';
D         = W' * exp(1j * PHI);
end

function L = kl_obj(Ry, R_hat, p_vec, lam_reg)
try
    L = real(log(det(Ry + 1e-12*eye(size(Ry,1)))) + trace(Ry \ R_hat)) ...
        + lam_reg * sum(p_vec);
catch
    L = Inf;
end
end
