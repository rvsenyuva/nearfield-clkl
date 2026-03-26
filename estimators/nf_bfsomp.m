function [theta_hat, r_hat, info] = nf_bfsomp(X, W, P)
%NF_BFSOMP  Beam-Focused Simultaneous OMP (BF-SOMP).
%
%  Implements Hussain, Abdallah & Eltawil,
%  "Redefining Polar Boundaries for Near-Field Channel Estimation
%  for Ultra-Massive MIMO Antenna Array,"
%  IEEE Trans. Wireless Commun., vol.24, no.10, Oct. 2025.
%  Algorithm 3 (BF-SOMP), adapted to our narrowband single-carrier model.
%
%  KEY DIFFERENCES FROM P-SOMP:
%  -----------------------------
%  1. POLAR CODEBOOK: same EBRD-bounded beam-depth sampling as nf_psomp v3
%     (Hussain Algorithm 1), but here the codebook Phi is the FULL ARRAY
%     steering matrix (M x S), not the compressed dictionary (N_RF x S).
%     BF-SOMP operates on Y = W^H*X (compressed snapshots), not R_hat.
%
%  2. ADAPTIVE SPARSITY: BF-SOMP does NOT require d to be known a priori.
%     It iterates until the spectral norm of the residual falls below:
%       epsilon_sigma = sigma * sqrt(M + 2*sqrt(M*log(M)))
%     where sigma = sqrt(N0) is the noise standard deviation.
%     P.d is used only as a safety upper cap.
%
%  3. WHITENING: noise covariance C = N0*(W^H*W) is Cholesky-factored and
%     applied to both the observations and the dictionary before SOMP,
%     ensuring the orthogonality conditions of OMP are met even when
%     W^H*W != I (as in our random-phase combiner).
%
%  4. CHANNEL RECONSTRUCTION: H_hat = Phi * H_rho (direct), where H_rho
%     is the sparse channel in the polar domain.  This differs from the
%     least-squares reconstruction used by CL-KL and P-SOMP.
%
%  INPUT/OUTPUT:
%  X         -- M x N full array snapshot matrix (same as DL-OMP, DFrFT)
%  W         -- M x N_RF hybrid combiner
%  P         -- parameter struct from nf_params()
%  theta_hat -- d_hat x 1 estimated angles [rad]
%  r_hat     -- d_hat x 1 estimated ranges [m]
%  info      -- diagnostic struct (see end of function)
%
%  NARROWBAND ADAPTATION NOTES:
%  The Hussain paper assumes OFDM with K subcarriers.  In our single-
%  carrier model we have N snapshots.  The multi-carrier sum in the
%  SOMP correlation step becomes a sum over N snapshots:
%    c_alpha = sum_{n=1}^{N} Psi_bar^H * r_bar(:,n)
%  which is equivalent to Hussain Algorithm 3 lines 12 with K=N.
%  The residual norm check uses ||R_bar||_2 (spectral norm of the
%  N_RF x N residual matrix).

M    = P.M;
N    = P.N;
N_RF = P.N_RF;
d    = P.d;        % used as UPPER CAP only (BF-SOMP estimates d adaptively)
lambda = P.lambda;
d_ant  = P.d_ant;

D_ap  = (M-1)*d_ant;
r_RD  = 2*D_ap^2/lambda;
r_min = P.r_lo_fac * r_RD;
r_max = P.r_hi_fac * r_RD;

Q_theta    = P.Q_theta;
theta_grid = linspace(P.theta_lo, P.theta_hi, Q_theta);
m_bar      = ((0:M-1).' - (M-1)/2);   % M x 1
m_bar2     = m_bar.^2;

% ====================================================================
%  STEP 1: Build EBRD-bounded beam-depth polar codebook Phi (M x S)
%          (Hussain Algorithm 1 -- same as nf_psomp v3 but uncompressed)
% ====================================================================
Phi_cell    = cell(1, Q_theta);
theta_atoms = [];
r_atoms     = [];

for ni = 1:Q_theta
    theta_n  = theta_grid(ni);
    cos2_n   = cos(theta_n)^2;
    r_EBRD_n = r_RD/10 * cos2_n;

    % r_BD-based range sampling (Hussain Theorem 1)
    r_samps = r_min;
    r_F     = r_min;
    n_s     = 1;
    while r_F < r_max && n_s < 8
        if r_F >= r_EBRD_n * 0.99; break; end
        denom = r_RD*cos2_n*(1 - 100*r_F^2/(r_RD^2*cos2_n^2));
        if denom <= 1e-12; break; end
        r_BD_n  = 20*r_F^2/denom;
        r_next  = min(r_F + r_BD_n, r_max);
        if abs(r_next - r_samps(end)) > 0.01
            r_samps(end+1) = r_next; %#ok<AGROW>
            n_s = n_s + 1;
        end
        if abs(r_next - r_max) < 1e-4; break; end
        r_F = r_next;
    end
    if abs(r_samps(end) - r_max) > 0.01
        r_samps(end+1) = r_max;
    end
    % Min 4 samples per angle
    if numel(r_samps) < 4
        n_add = 4 - numel(r_samps);
        for ai = 1:n_add
            extra = r_min + (r_max-r_min)*ai/(n_add+1);
            if ~any(abs(extra - r_samps) < 0.1)
                r_samps(end+1) = extra; %#ok<AGROW>
            end
        end
        r_samps = sort(r_samps);
    end

    Qr_n    = numel(r_samps);
    omega_n = (2*pi*d_ant/lambda)*cos(theta_n);
    c_n     = (pi*d_ant^2/lambda)*sin(theta_n)^2;

    A_n = zeros(M, Qr_n);
    for ri = 1:Qr_n
        kappa = c_n / r_samps(ri);
        A_n(:,ri) = exp(1j*omega_n*m_bar - 1j*kappa*m_bar2);
    end
    Phi_cell{ni}  = A_n;   % M x Qr_n  (uncompressed)
    theta_atoms   = [theta_atoms, repmat(theta_n, 1, Qr_n)]; %#ok<AGROW>
    r_atoms       = [r_atoms,     r_samps]; %#ok<AGROW>
end

Phi    = cat(2, Phi_cell{:});   % M x S  (full-array polar codebook)
S      = size(Phi, 2);

% ====================================================================
%  STEP 2: Compressed observations and noise threshold
% ====================================================================
Y_comp = W' * X;   % N_RF x N

% Noise power estimate via eigenvalue decomposition of sample covariance
R_x    = (1/N) * (X * X');
ev     = sort(real(eig((R_x+R_x')/2)), 'ascend');
n_noise_ev = max(1, M - d);
N0_est = max(mean(ev(1:n_noise_ev)), 1e-12);
sigma  = sqrt(N0_est);

% Stopping threshold (Hussain eq., line 4 of Algorithm 3)
eps_sigma = sigma * sqrt(M + 2*sqrt(M*log(max(M,2))));

% ====================================================================
%  STEP 3: Whitening (Hussain Algorithm 3, lines 1-3)
%  Noise covariance: C = N0 * (W^H * W)
%  Cholesky: C = D_ch * D_ch^H
%  Whitened observations: Y_bar   = D_ch^{-1} * Y_comp
%  Whitened sensing:      Psi_bar = D_ch^{-1} * W^H * Phi
% ====================================================================
C_noise = N0_est * (W' * W);
C_noise = (C_noise + C_noise')/2 + 1e-10*eye(N_RF);   % symmetrise + regularise
try
    D_ch = chol(C_noise, 'lower');    % lower-triangular: C = D*D^H
    D_ch_inv = inv(D_ch);
catch
    D_ch_inv = eye(N_RF) / max(sqrt(N0_est), 1e-8);
end

Y_bar   = D_ch_inv * Y_comp;          % N_RF x N  (whitened observations)
WPhi    = W' * Phi;                    % N_RF x S  (compressed codebook)
Psi_bar = D_ch_inv * WPhi;            % N_RF x S  (whitened sensing matrix)

% ====================================================================
%  STEP 4: Adaptive SOMP (Hussain Algorithm 3, lines 8-20)
% ====================================================================
R_bar   = Y_bar;       % N_RF x N  (residual, initialised to Y_bar)
gamma   = [];          % support set
d_hat   = 0;

while norm(R_bar, 2) > eps_sigma && d_hat < d
    % Correlation across all N snapshots (line 12)
    c_alpha = sum(Psi_bar' * R_bar, 2);     % S x 1  (real part drives selection)
    [~, alpha_star] = max(abs(c_alpha));

    % Update support
    gamma(end+1) = alpha_star; %#ok<AGROW>
    d_hat = d_hat + 1;

    % Orthogonal projection onto support (line 15)
    Psi_gamma = Psi_bar(:, gamma);          % N_RF x d_hat
    H_rho_gamma = pinv(Psi_gamma) * Y_bar;  % d_hat x N

    % Update residual (line 16)
    R_bar = Y_bar - Psi_gamma * H_rho_gamma;
end

% Fallback: if nothing selected, take top d atoms
if isempty(gamma)
    c_alpha = sum(Psi_bar' * R_bar, 2);
    [~, idx] = sort(abs(c_alpha), 'descend');
    gamma = idx(1:min(d, S)).';
    d_hat = numel(gamma);
end

% ====================================================================
%  STEP 5: Extract angle and range estimates
% ====================================================================
theta_hat = theta_atoms(gamma).';
r_hat     = r_atoms(gamma).';

% Pad to d if fewer paths found
if numel(theta_hat) < d
    n_miss    = d - numel(theta_hat);
    theta_hat = [theta_hat; mean([P.theta_lo, P.theta_hi])*ones(n_miss,1)];
    r_hat     = [r_hat;     0.5*(r_min+r_max)*ones(n_miss,1)];
end
theta_hat = theta_hat(1:d);
r_hat     = r_hat(1:d);

% ====================================================================
info.d_hat        = d_hat;
info.gamma        = gamma;
info.eps_sigma    = eps_sigma;
info.N0_est       = N0_est;
info.Q_total      = S;
info.sampling     = 'beam_depth_rBD_BF-SOMP';
info.norm_R_final = norm(R_bar, 2);
end
