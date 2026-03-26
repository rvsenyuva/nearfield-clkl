function [rmse_theta_crb, rmse_r_crb] = nf_crb(theta_true, r_true, p_true, N0, W, P)
%NF_CRB  Stochastic Cramer-Rao Bound for compressed near-field observations.
%
%  Computes the CRB for angle (theta) and range (r) estimation from the
%  COMPRESSED observation model:
%
%    y_n = W^H * (A*s_n + w_n),   n = 1..N
%    s_n ~ CN(0, diag(p)),  w_n ~ CN(0, N0*I_M)
%    y_n ~ CN(0, R_y),   R_y = W^H*(A*diag(p)*A^H + N0*I_M)*W
%
%  where W in C^{M x N_RF} is the hybrid combiner and A in C^{M x d} is
%  the Fresnel steering matrix.
%
%  FORMULA (Stoica & Nehorai 1990, eq. 3.4; Ottersten et al. 1998):
%    FIM_{ij} = N * Re( tr(R_y^{-1} dR_y/deta_i  R_y^{-1} dR_y/deta_j) )
%
%  PARAMETER VECTOR:
%    eta = [omega_1..omega_d, kappa_1..kappa_d, p_1..p_d, N0]   (3d+1 params)
%    omega_l = (2*pi*d/lambda)*cos(theta_l)          (spatial frequency)
%    kappa_l = (pi*d^2/lambda)*sin^2(theta_l) / r_l  (Fresnel curvature)
%
%  DERIVATIVES (only non-zero columns of dA shown):
%    dA/d_omega_l  : column l =  j * m_bar .* a_l,    rest zero
%    dA/d_kappa_l  : column l = -j * m_bar^2 .* a_l,  rest zero
%    dR_y/d_omega_l = W^H*(dA/do_l*P*A^H + A*P*dA^H/do_l)*W
%    dR_y/d_kappa_l = W^H*(dA/dk_l*P*A^H + A*P*dA^H/dk_l)*W
%    dR_y/d_p_l     = W^H * a_l * a_l^H * W
%    dR_y/d_N0      = W^H * W
%
%  ERROR PROPAGATION to physical parameters:
%    d_omega/d_theta_l = -(2*pi*d/lambda)*sin(theta_l)
%    d_kappa/d_r_l     = -(pi*d^2/lambda)*sin^2(theta_l) / r_l^2
%    CRB_theta_l = [J^{-1}]_{ll}         / (d_omega/d_theta_l)^2
%    CRB_r_l     = [J^{-1}]_{d+l, d+l}  / (d_kappa/d_r_l)^2
%
%  INPUTS:
%    theta_true  -- d x 1 true angles [rad]
%    r_true      -- d x 1 true ranges [m]
%    p_true      -- d x 1 true path powers (normalised)
%    N0          -- true noise variance (scalar)
%    W           -- M x N_RF hybrid combiner
%    P           -- parameter struct from nf_params()
%
%  OUTPUTS:
%    rmse_theta_crb -- sqrt of average CRB across d paths [deg]
%    rmse_r_crb     -- sqrt of average CRB across d paths [m]
%
%  NOTE: This is the COMPRESSED-DOMAIN CRB (N_RF x N_RF covariance).
%  It accounts for information loss due to hybrid combining and is strictly
%  larger than the full-array CRB.  The full-array CRB (Grosicki et al.
%  2005) can be obtained by setting W = I_M, but is only relevant for
%  fully-digital architectures.
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

M      = P.M;
N      = P.N;
N_RF   = size(W,2);
d      = numel(theta_true);
lambda = P.lambda;
d_ant  = P.d_ant;
n_par  = 3*d + 1;   % omega(1..d), kappa(1..d), p(1..d), N0

m_bar  = ((0:M-1).' - (M-1)/2);   % M x 1
m_bar2 = m_bar.^2;

% ---- Build steering matrix A and compressed covariance R_y -----------
A = zeros(M, d);
omega_vec = zeros(d,1);
kappa_vec = zeros(d,1);
c_coef    = zeros(d,1);  % (pi*d_ant^2/lambda)*sin^2(theta_l)

for ell = 1:d
    omega_vec(ell) = (2*pi*d_ant/lambda)*cos(theta_true(ell));
    c_coef(ell)    = (pi*d_ant^2/lambda)*sin(theta_true(ell))^2;
    kappa_vec(ell) = c_coef(ell)/r_true(ell);
    lin_ph  = omega_vec(ell) * m_bar;
    quad_ph = kappa_vec(ell) * m_bar2;
    A(:,ell) = exp(1j*lin_ph - 1j*quad_ph);
end

WA    = W' * A;                             % N_RF x d
R_h   = WA * diag(p_true) * WA';           % N_RF x N_RF  (signal term)
R_n   = N0 * (W'*W);                       % N_RF x N_RF  (noise term)
R_y   = R_h + R_n;
R_y   = (R_y + R_y')/2 + 1e-10*eye(N_RF); % symmetrise + small regularisation
R_inv = R_y \ eye(N_RF);                   % N_RF x N_RF (preferred over inv)

% ---- Build all 3d+1 derivative matrices dR_y/d_eta_i -----------------
dRy = cell(n_par, 1);

for ell = 1:d
    a_ell  = A(:,ell);
    Wa_ell = WA(:,ell);    % W' * a_ell

    % dR_y/d_omega_l
    da_do = 1j * m_bar .* a_ell;              % M x 1
    Wda   = W' * da_do;                        % N_RF x 1
    dRy{ell} = p_true(ell) * ...
        (Wda*Wa_ell' + Wa_ell*Wda');

    % dR_y/d_kappa_l
    da_dk  = -1j * m_bar2 .* a_ell;           % M x 1
    Wdak   = W' * da_dk;                       % N_RF x 1
    dRy{d+ell} = p_true(ell) * ...
        (Wdak*Wa_ell' + Wa_ell*Wdak');

    % dR_y/d_p_l
    dRy{2*d+ell} = Wa_ell * Wa_ell';
end

% dR_y/d_N0
dRy{3*d+1} = W' * W;

% ---- Fisher Information Matrix FIM (3d+1 x 3d+1) --------------------
%  FIM_{ij} = N * Re( tr(R_inv * dRy_i * R_inv * dRy_j) )
FIM = zeros(n_par, n_par);
for ii = 1:n_par
    RiDi = R_inv * dRy{ii};      % precompute once per row
    for jj = ii:n_par             % exploit symmetry
        FIM(ii,jj) = N * real(trace(RiDi * R_inv * dRy{jj}));
        FIM(jj,ii) = FIM(ii,jj);
    end
end

% ---- Invert FIM via SVD-based pseudoinverse (robust to near-singularity)
%  Tikhonov \-solve can still blow up for rank-deficient FIMs.
%  SVD pinv with relative tolerance gracefully zeroes out near-zero
%  singular values, giving a finite (though loose) CRB for all cases.
%  Tolerance: 1e-6 * sigma_max  (keeps singular values down to 1e-6
%  of the largest, which corresponds to ~10^6 condition number).
[U_fim, S_fim, V_fim] = svd(FIM);
sv      = diag(S_fim);
tol_sv  = 1e-6 * sv(1);                     % relative tolerance
sv_inv  = zeros(size(sv));
sv_inv(sv > tol_sv) = 1 ./ sv(sv > tol_sv); % invert only above tol
CRB_mat = V_fim * diag(sv_inv) * U_fim';     % CRB = FIM^{-1} (pseudo)

% ---- Extract marginal CRBs and apply error propagation ---------------
crb_theta = zeros(d,1);
crb_r     = zeros(d,1);

for ell = 1:d
    % CRB for omega_ell (index ell in parameter vector)
    var_omega = max(CRB_mat(ell, ell), 0);

    % CRB for kappa_ell (index d+ell)
    var_kappa = max(CRB_mat(d+ell, d+ell), 0);

    % Error propagation: omega -> theta
    %   d_omega/d_theta = -(2*pi*d_ant/lambda)*sin(theta)
    domega_dtheta = -(2*pi*d_ant/lambda) * sin(theta_true(ell));
    crb_theta(ell) = var_omega / domega_dtheta^2;

    % Error propagation: kappa -> r
    %   d_kappa/d_r = -c_coef(ell)/r^2
    dkappa_dr     = -c_coef(ell) / r_true(ell)^2;
    crb_r(ell)    = var_kappa / dkappa_dr^2;
end

% ---- Average RMSE over paths -----------------------------------------
rmse_theta_crb = mean(sqrt(crb_theta)) * 180/pi;   % convert to degrees
rmse_r_crb     = mean(sqrt(crb_r));                 % in metres
end
