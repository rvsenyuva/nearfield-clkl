function [out1, out2] = nf_crb(theta_true, r_true, p_true, N0, W, P, mode)
%NF_CRB  Stochastic Cramér-Rao Bound for near-field observations.
%
%  Computes the CRB for angle (theta) and range (r) estimation under
%  one of two observation models, selected by the optional MODE argument.
%
%  ---- COMPRESSED-DOMAIN MODEL (mode = 'compressed', DEFAULT) ----------
%
%    y_n = W^H*(A*s_n + w_n),   n = 1..N
%    y_n ~ CN(0, R_y),   R_y = W^H*(A*diag(p)*A^H + N0*I_M)*W
%
%  Correct bound for methods observing only R_hat_y (CL-KL, P-SOMP).
%
%  ---- FULL-ARRAY MODEL (mode = 'full') --------------------------------
%
%    x_n = A*s_n + w_n,   n = 1..N,   x_n ~ CN(0, R_x)
%    R_x = A*diag(p)*A^H + N0*I_M
%
%  Obtained by substituting W = I_M in the compressed model.  Correct
%  bound for fully-digital receivers (DL-OMP, MUSIC+Tri, DFrFT-NOMP,
%  BF-SOMP).  Data-processing inequality: CRB_full <= CRB_comp, with
%  strict inequality whenever N_RF < M.
%
%  ---- BOTH MODES (mode = 'both') --------------------------------------
%
%  Returns a struct in out1 with fields:
%    out1.theta_comp  -- compressed CRB for angle [deg, sqrt of avg CRB]
%    out1.r_comp      -- compressed CRB for range [m,   sqrt of avg CRB]
%    out1.theta_full  -- full-array CRB for angle [deg, sqrt of avg CRB]
%    out1.r_full      -- full-array CRB for range [m,   sqrt of avg CRB]
%  out2 is unused in 'both' mode.
%
%  ---- SCALAR CALL SYNTAX (mode = 'compressed' or 'full') --------------
%
%    [crb_theta_deg, crb_r_m] = nf_crb(theta, r, p, N0, W, P)
%    [crb_theta_deg, crb_r_m] = nf_crb(theta, r, p, N0, W, P, 'compressed')
%    [crb_theta_deg, crb_r_m] = nf_crb(theta, r, p, N0, W, P, 'full')
%
%  The default ('compressed') call is bit-identical to all pre-P9
%  callers: nargout, argument count, and return values are unchanged.
%
%  ---- FORMULA (Stoica & Nehorai 1990, Eq. 3.4) -----------------------
%
%    FIM_{ij} = N * Re( tr(R^{-1} dR/deta_i  R^{-1} dR/deta_j) )
%
%    eta = [omega_1..omega_d, kappa_1..kappa_d, p_1..p_d, N0]  (3d+1)
%    omega_l = (2*pi*d_ant/lambda)*cos(theta_l)
%    kappa_l = (pi*d_ant^2/lambda)*sin^2(theta_l) / r_l
%
%  DERIVATIVES (compressed; full-array drops all W^H prefactors):
%    dR_y/d_omega_l = W^H*(dA/do_l*P*A^H + A*P*dA^H/do_l)*W
%    dR_y/d_kappa_l = W^H*(dA/dk_l*P*A^H + A*P*dA^H/dk_l)*W
%    dR_y/d_p_l     = W^H * a_l * a_l^H * W
%    dR_y/d_N0      = W^H * W
%
%  ERROR PROPAGATION to physical parameters:
%    d_omega/d_theta_l = -(2*pi*d_ant/lambda)*sin(theta_l)
%    d_kappa/d_r_l     = -(pi*d_ant^2/lambda)*sin^2(theta_l) / r_l^2
%    CRB_theta_l = [FIM^{-1}]_{ll}          / (d_omega/d_theta_l)^2
%    CRB_r_l     = [FIM^{-1}]_{d+l, d+l}   / (d_kappa/d_r_l)^2
%
%  INPUTS:
%    theta_true  -- d x 1 true angles [rad]
%    r_true      -- d x 1 true ranges [m]
%    p_true      -- d x 1 true path powers (normalised)
%    N0          -- true noise variance (scalar)
%    W           -- M x N_RF hybrid combiner (used only in 'compressed')
%    P           -- parameter struct from nf_params()
%    mode        -- 'compressed' (default) | 'full' | 'both'
%
%  OUTPUTS (mode = 'compressed' or 'full'):
%    out1  -- sqrt(mean CRB_theta) across d paths [deg]
%    out2  -- sqrt(mean CRB_r)     across d paths [m]
%
%  OUTPUT (mode = 'both'):
%    out1  -- struct with fields theta_comp, r_comp, theta_full, r_full
%    out2  -- [] (unused)
%
%  REFERENCES:
%  [1] P. Stoica & A. Nehorai, "Performance study of conditional and
%      unconditional direction-of-arrival estimation," IEEE Trans. ASSP,
%      vol.38, no.10, pp.1783-1795, Oct. 1990.
%  [2] B. Ottersten, P. Stoica & R. Roy, "Covariance matching estimation
%      techniques for array signal processing applications," Digit. Signal
%      Process., vol.8, no.3, pp.185-210, 1998.
%  [3] E. Grosicki, K. Abed-Meraim & Y. Hua, "A weighted linear prediction
%      method for near-field source localization," IEEE Trans. Signal
%      Process., vol.53, no.10, pp.3651-3660, Oct. 2005.

% ---- Parse mode argument (default: 'compressed') ---------------------
if nargin < 7 || isempty(mode)
    mode = 'compressed';
end
mode = lower(strtrim(mode));
if ~any(strcmp(mode, {'compressed','full','both'}))
    error('nf_crb: mode must be ''compressed'', ''full'', or ''both''. Got: ''%s''.', mode);
end

% ---- Shared parameters -----------------------------------------------
M      = P.M;
N      = P.N;
d      = numel(theta_true);
lambda = P.lambda;
d_ant  = P.d_ant;
n_par  = 3*d + 1;   % [omega(1..d), kappa(1..d), p(1..d), N0]

m_bar  = ((0:M-1).' - (M-1)/2);   % M x 1, centred element indices
m_bar2 = m_bar.^2;

% ---- Build Fresnel steering matrix A (M x d) -------------------------
A      = zeros(M, d);
c_coef = zeros(d,1);   % (pi*d_ant^2/lambda)*sin^2(theta_l) for error prop

for ell = 1:d
    omega_ell    = (2*pi*d_ant/lambda)*cos(theta_true(ell));
    c_coef(ell)  = (pi*d_ant^2/lambda)*sin(theta_true(ell))^2;
    kappa_ell    = c_coef(ell) / r_true(ell);
    A(:,ell)     = exp(1j*omega_ell*m_bar - 1j*kappa_ell*m_bar2);
end

% ---- Compute CRB(s) --------------------------------------------------
switch mode

    case 'compressed'
        [crb_th, crb_r] = branch_compressed( ...
            A, W, d, n_par, N, m_bar, m_bar2, p_true, N0, ...
            d_ant, lambda, theta_true, c_coef, r_true);
        out1 = mean(sqrt(crb_th)) * 180/pi;   % deg
        out2 = mean(sqrt(crb_r));              % m

    case 'full'
        [crb_th, crb_r] = branch_fullarray( ...
            A, d, n_par, N, m_bar, m_bar2, p_true, N0, M, ...
            d_ant, lambda, theta_true, c_coef, r_true);
        out1 = mean(sqrt(crb_th)) * 180/pi;
        out2 = mean(sqrt(crb_r));

    case 'both'
        [crb_th_c, crb_r_c] = branch_compressed( ...
            A, W, d, n_par, N, m_bar, m_bar2, p_true, N0, ...
            d_ant, lambda, theta_true, c_coef, r_true);
        [crb_th_f, crb_r_f] = branch_fullarray( ...
            A, d, n_par, N, m_bar, m_bar2, p_true, N0, M, ...
            d_ant, lambda, theta_true, c_coef, r_true);
        s.theta_comp = mean(sqrt(crb_th_c)) * 180/pi;
        s.r_comp     = mean(sqrt(crb_r_c));
        s.theta_full = mean(sqrt(crb_th_f)) * 180/pi;
        s.r_full     = mean(sqrt(crb_r_f));
        out1 = s;
        out2 = [];
end
end   % nf_crb


% =====================================================================
%  BRANCH: compressed-domain CRB
%  R_y = W^H*(A*diag(p)*A^H + N0*I_M)*W   (N_RF x N_RF)
% =====================================================================
function [crb_theta, crb_r] = branch_compressed( ...
        A, W, d, n_par, N, m_bar, m_bar2, p_true, N0, ...
        d_ant, lambda, theta_true, c_coef, r_true)

N_RF  = size(W,2);
WA    = W' * A;                              % N_RF x d

R_h   = WA * diag(p_true) * WA';
R_n   = N0 * (W'*W);
R_y   = R_h + R_n;
R_y   = (R_y + R_y')/2 + 1e-10*eye(N_RF);  % symmetrise + regularise
R_inv = R_y \ eye(N_RF);

% Build derivative matrices dR_y/d_eta_i  (N_RF x N_RF each)
dR = cell(n_par, 1);
for ell = 1:d
    a_ell  = A(:,ell);
    Wa_ell = WA(:,ell);

    da_do = 1j * m_bar .* a_ell;
    Wda   = W' * da_do;
    dR{ell} = p_true(ell) * (Wda*Wa_ell' + Wa_ell*Wda');

    da_dk = -1j * m_bar2 .* a_ell;
    Wdak  = W' * da_dk;
    dR{d+ell} = p_true(ell) * (Wdak*Wa_ell' + Wa_ell*Wdak');

    dR{2*d+ell} = Wa_ell * Wa_ell';         % dR_y/d_p_l
end
dR{3*d+1} = W' * W;                         % dR_y/d_N0

[crb_theta, crb_r] = fim_to_crb(dR, R_inv, N, n_par, d, ...
    d_ant, lambda, theta_true, c_coef, r_true);
end


% =====================================================================
%  BRANCH: full-array CRB
%  R_x = A*diag(p)*A^H + N0*I_M   (M x M)
%  Specialisation W = I_M: W^H*(..)*W -> (..), W^H*a_l -> a_l,
%                           W^H*W -> I_M.
% =====================================================================
function [crb_theta, crb_r] = branch_fullarray( ...
        A, d, n_par, N, m_bar, m_bar2, p_true, N0, M, ...
        d_ant, lambda, theta_true, c_coef, r_true)

R_x   = A * diag(p_true) * A' + N0 * eye(M);
R_x   = (R_x + R_x')/2 + 1e-10*eye(M);     % symmetrise + regularise
R_inv = R_x \ eye(M);

% Build derivative matrices dR_x/d_eta_i  (M x M each)
% W^H prefactors drop out: da_do stays M x 1, no W' projection.
dR = cell(n_par, 1);
for ell = 1:d
    a_ell = A(:,ell);

    da_do = 1j * m_bar .* a_ell;
    dR{ell} = p_true(ell) * (da_do*a_ell' + a_ell*da_do');

    da_dk = -1j * m_bar2 .* a_ell;
    dR{d+ell} = p_true(ell) * (da_dk*a_ell' + a_ell*da_dk');

    dR{2*d+ell} = a_ell * a_ell';           % dR_x/d_p_l
end
dR{3*d+1} = eye(M);                         % dR_x/d_N0 (W^H*W = I_M)

[crb_theta, crb_r] = fim_to_crb(dR, R_inv, N, n_par, d, ...
    d_ant, lambda, theta_true, c_coef, r_true);
end


% =====================================================================
%  SHARED: FIM assembly, SVD pseudoinverse, error propagation
% =====================================================================
function [crb_theta, crb_r] = fim_to_crb(dR, R_inv, N, n_par, d, ...
        d_ant, lambda, theta_true, c_coef, r_true)
%FIM_TO_CRB  Assemble FIM, invert, extract marginal CRBs.

% FIM_{ij} = N * Re( tr(R^{-1} dR_i  R^{-1} dR_j) )
FIM = zeros(n_par, n_par);
for ii = 1:n_par
    RiDi = R_inv * dR{ii};           % precompute R^{-1}*dR_i
    for jj = ii:n_par                 % upper triangle only (exploit symmetry)
        FIM(ii,jj) = N * real(trace(RiDi * R_inv * dR{jj}));
        FIM(jj,ii) = FIM(ii,jj);
    end
end

% SVD pseudoinverse with relative threshold 1e-6*sigma_max.
% Gracefully handles rank-deficient FIMs near the identifiability boundary
% (d > floor((N_RF-1)/2) for compressed model) without producing Inf/NaN.
[U_fim, S_fim, V_fim] = svd(FIM);
sv      = diag(S_fim);
tol_sv  = 1e-6 * sv(1);
sv_inv  = zeros(size(sv));
sv_inv(sv > tol_sv) = 1 ./ sv(sv > tol_sv);
CRB_mat = V_fim * diag(sv_inv) * U_fim';

% Marginal CRBs with error propagation
crb_theta = zeros(d,1);
crb_r     = zeros(d,1);
for ell = 1:d
    var_omega = max(CRB_mat(ell,   ell),   0);
    var_kappa = max(CRB_mat(d+ell, d+ell), 0);

    domega_dtheta  = -(2*pi*d_ant/lambda) * sin(theta_true(ell));
    crb_theta(ell) = var_omega / domega_dtheta^2;

    dkappa_dr   = -c_coef(ell) / r_true(ell)^2;
    crb_r(ell)  = var_kappa / dkappa_dr^2;
end
end
