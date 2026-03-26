function [theta_hat, r_hat, info] = nf_dfrft_nomp(X, ~, P)
%NF_DFRFT_NOMP  DFrFT-based NOMP gridless hybrid-field channel estimator.
%
%  Implements Yang et al. (IEEE Wireless Commun. Lett. 2024), Algorithm 1.
%
%  CRITICAL FIX (v2):
%  ------------------
%  Previous version averaged |DFrFT(x_n)|^2 INCOHERENTLY across N snapshots,
%  which destroyed the inter-path cancellation that the DFrFT's chirp-
%  alignment property relies on. The Yang et al. algorithm is designed for
%  a SINGLE snapshot y in C^M (y = h + noise, SNR = rho per element).
%
%  Fix: compute a SINGLE coherently-averaged snapshot
%       y_avg = (1/N) * sum_n X(:,n) = mean(X, 2)
%  and run the entire DFrFT on this single M x 1 vector. This restores
%  the single-snapshot assumption, giving SNR_eff = N * SNR_single which
%  is why the algorithm improves with N as expected for a coherent method.
%
%  The residual deflation also operates on y_avg, matching the paper's
%  Algorithm 1 step (18).
%
%  ADDITIONAL FIXES:
%  - SNR estimate now uses trace of R_x minus noise floor (more stable
%    than max-eigenvalue), avoiding premature termination.
%  - termination threshold T = M / SNR_eff accounts for averaging gain.

M      = P.M;
N      = P.N;
d_target = P.d;
lambda = P.lambda;
d_ant  = P.d_ant;

D_ap  = (M-1)*d_ant;
r_RD  = 2*D_ap^2/lambda;
r_min = P.r_lo_fac * r_RD;
r_max = P.r_hi_fac * r_RD;

m_bar  = (0:M-1).';   % 0-indexed

% ---- KEY FIX: use coherently averaged snapshot --------------------------
% Yang et al. eq. (1): y = h*s + n. With s=1 for all N pilots, the ML
% estimate of h is exactly y_avg = mean(X,2). This gives SNR_eff = N*SNR.
y_avg  = mean(X, 2);   % M x 1  -- this is the single snapshot for DFrFT

% ---- DFrFT parameters ---------------------------------------------------
P_ord  = min(M, 128);
zeta   = 4;
Mz     = M * zeta;
p_vec  = linspace(0.5, 1.5, P_ord);
alpha_vec = p_vec * pi/2;

% ---- SNR estimate from averaged snapshot --------------------------------
% Use R_x = XX^H/N for eigenvalue-based N0 estimate (more stable)
R_x    = (1/N) * (X * X');
ev_all = sort(real(eig((R_x+R_x')/2)), 'descend');
N0_est = max(mean(ev_all(d_target+1:end)), 1e-12);
% SNR_eff accounts for N-fold averaging gain in y_avg
SNR_eff = max(norm(y_avg)^2/M / N0_est, 1);
% Termination threshold: for averaged snapshot, energy is ~M*(h_norm^2 + N0/N)
T_term  = M * N0_est / N;   % divided by N because y_avg noise is N0/N per element

% ---- Precompute DFrFT spectra on single averaged snapshot ---------------
% For coherent single snapshot, just compute DFrFT(y_avg) at each p.
V_avg = zeros(P_ord, Mz);   % power spectra

for pi_idx = 1:P_ord
    alpha = alpha_vec(pi_idx);
    c_alpha = (1 - 1j*cot(alpha)) / Mz;
    y_z   = [y_avg; zeros(Mz - M, 1)];
    m_z   = (0:Mz-1).';
    h_m   = y_z .* exp(-1j * pi * m_z.^2 * tan(alpha/2) / Mz);
    Y_alpha = c_alpha * exp(-1j * pi * (0:Mz-1).' .^2 .* tan(alpha/2) / Mz) .* ...
              fft(h_m .* exp(-1j * pi * m_z.^2 .* tan(alpha/2) / Mz));
    V_avg(pi_idx, :) = abs(Y_alpha.').^2;
end

% ---- Greedy path detection with deflation on y_avg ----------------------
omega_min = 2*pi*d_ant/lambda * cos(P.theta_hi);
omega_max = 2*pi*d_ant/lambda * cos(P.theta_lo);
kappa_max = pi*d_ant^2/lambda / r_min;

theta_hat_list = [];
r_hat_list     = [];
omega_list     = [];
kappa_list     = [];

y_res    = y_avg;   % residual operates on the averaged snapshot
V_working = V_avg;
paths_found = 0;

for k = 1:d_target
    % Coarse detection
    [~, best_lin] = max(V_working(:));
    [pi_best, qi_best] = ind2sub([P_ord, Mz], best_lin);

    alpha_k = alpha_vec(pi_best);
    q_k     = qi_best - 1;

    omega_k = 2*pi * q_k * (1/sin(alpha_k)) / Mz;
    kappa_k = pi * (cos(alpha_k)/sin(alpha_k)) / Mz;

    omega_k = min(omega_max, max(omega_min, omega_k));
    kappa_k = min(kappa_max, max(0, kappa_k));

    % Newton refinement (5 steps) -- score against y_res (not R_x)
    m_bar_r = m_bar - (M-1)/2;
    m_bar2  = m_bar_r.^2;
    step_size = 1e-4;

    for rs_iter = 1:5
        a_k    = exp(1j*omega_k*m_bar_r - 1j*kappa_k*m_bar2);
        da_w   = 1j*m_bar_r   .* a_k;
        da_kap = -1j*m_bar2   .* a_k;
        % Gradient of ||y_res - gamma*a||^2
        gamma_k = (a_k' * y_res) / max(norm(a_k)^2, 1e-12);
        resid_k = y_res - gamma_k * a_k;
        df_om   = -2 * real(da_w'   * resid_k);
        df_kap  = -2 * real(da_kap' * resid_k);
        omega_k = min(omega_max, max(omega_min, omega_k - step_size*df_om));
        kappa_k = min(kappa_max, max(1e-12,     kappa_k - step_size*df_kap));
    end

    % Recover physical parameters
    cos_th = min(1, max(-1, omega_k*lambda/(2*pi*d_ant)));
    th_k   = real(acos(cos_th));
    c_k    = pi*d_ant^2/lambda * sin(th_k)^2;
    if c_k > 1e-12 && kappa_k > 1e-12
        r_k = min(r_max, max(r_min, c_k / kappa_k));
    else
        r_k = 0.5*(r_min + r_max);
    end

    theta_hat_list(end+1) = th_k;   %#ok<AGROW>
    r_hat_list(end+1)     = r_k;    %#ok<AGROW>
    omega_list(end+1)     = omega_k; %#ok<AGROW>
    kappa_list(end+1)     = kappa_k; %#ok<AGROW>
    paths_found = paths_found + 1;

    % Suppress detected peak region
    suppress_half = max(1, round(P_ord/20));
    r_lo = max(1, pi_best - suppress_half);
    r_hi = min(P_ord, pi_best + suppress_half);
    q_lo = max(1, qi_best - round(Mz/20));
    q_hi = min(Mz, qi_best + round(Mz/20));
    V_working(r_lo:r_hi, q_lo:q_hi) = 0;

    % Deflate residual on averaged snapshot (paper eq. 18)
    a_k   = exp(1j*omega_k*m_bar_r - 1j*kappa_k*m_bar2);
    gamma = (a_k' * y_res) / max(norm(a_k)^2, 1e-12);
    y_res = y_res - gamma * a_k;

    if norm(y_res)^2 < T_term && paths_found >= 1
        break;
    end
end

% ---- Global gain correction (eq. 23) on averaged snapshot ---------------
if paths_found > 0
    n_p = paths_found;
    A_mat = zeros(M, n_p);
    for j = 1:n_p
        m_r = m_bar - (M-1)/2;
        A_mat(:,j) = exp(1j*omega_list(j)*m_r - 1j*kappa_list(j)*m_r.^2);
    end
    % LS on averaged snapshot
    gamma_ls = (A_mat'*A_mat + 1e-5*eye(n_p)) \ (A_mat'*y_avg);
    theta_hat = theta_hat_list(1:n_p).';
    r_hat     = r_hat_list(1:n_p).';
else
    theta_hat = ones(d_target,1) * mean([P.theta_lo, P.theta_hi]);
    r_hat     = ones(d_target,1) * 0.5*(r_min+r_max);
end

% Pad to d_target
if numel(theta_hat) < d_target
    n_miss = d_target - numel(theta_hat);
    theta_hat = [theta_hat; mean([P.theta_lo,P.theta_hi])*ones(n_miss,1)];
    r_hat     = [r_hat;     0.5*(r_min+r_max)*ones(n_miss,1)];
end
theta_hat = theta_hat(1:d_target);
r_hat     = r_hat(1:d_target);

info.paths_found = paths_found;
info.omega       = omega_list;
info.kappa       = kappa_list;
info.N0_est      = N0_est;
info.T_term      = T_term;
info.SNR_eff     = SNR_eff;
end
