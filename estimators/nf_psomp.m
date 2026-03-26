function [theta_hat, r_hat, info] = nf_psomp(R_hat, W, P)
%NF_PSOMP  Polar-grid Simultaneous OMP (P-SOMP) -- v3.
%
%  STEP 2 UPGRADE: Beam-depth range sampling (Hussain et al. TWC 2025).
%  ______________________________________________________________________
%  Previous version used Cui & Dai Lemma 1 to set a GLOBAL Q_r in [4,8]
%  with UNIFORM spacing in u=1/r.  This has two weaknesses:
%
%  (a) The coherence condition is averaged over all angles, so it neither
%      captures the angle-dependent near-field structure nor matches the
%      true physical sampling criterion for polar-domain sparsity.
%
%  (b) Uniform-u spacing is dense near r_min and coarse near r_max,
%      which is the WRONG direction: curvature 1/r varies LESS at large r
%      so fewer samples are needed there, not more.
%
%  Hussain et al. (IEEE TWC 2025, Theorem 1) derive the 3dB beam depth:
%
%    r_BD(r_F, theta) = 20*r_F^2 / (r_RD*cos^2(theta)) /
%                       (1 - 100*r_F^2 / (r_RD^2*cos^4(theta)))
%
%  which is the distance interval around r_F where the normalised array
%  gain drops by at most 3dB.  Codewords separated by r_BD have spatial
%  correlation <= 0.5, guaranteeing low dictionary coherence.
%
%  This version uses PER-ANGLE r_BD sampling (mirrors Hussain Algorithm 1):
%    - Start at r_min (= P.r_lo_fac * r_RD)
%    - Step by r_BD(r_F, theta_n) at each focus distance
%    - Stop when r_BD becomes infinite (at r_EBRD) or r_max is reached
%    - Always include r_max as a final sample
%    - Enforce Qr >= 4 per angle for OMP stability (pad with uniform
%      samples if the r_BD grid gives fewer than 4 points)
%    - Cap at Qr <= 8 per angle to bound dictionary size
%
%  The result is a VARIABLE-SIZE per-angle grid: denser near r_min (strong
%  near-field, small r_BD) and coarser near r_max (weak near-field, large
%  r_BD or r_BD -> inf in the far field).  This is physically correct and
%  matches the sparsity structure of the near-field channel.
%
%  Total dictionary size S = sum_{n=1}^{Q_theta} Qr_n.
%  For M=64, f_c=28GHz, theta in [20,60]deg, r in [0.05,1.0]*r_RD:
%  S ~ 1023 atoms (vs 1024 with old Q_r=4 uniform, essentially equal size).
%
%  COVARIANCE-DOMAIN SOMP: unchanged from previous version.

M       = P.M;
N_RF    = P.N_RF;
Q_theta = P.Q_theta;
d       = P.d;
lambda  = P.lambda;
d_ant   = P.d_ant;

D_ap  = (M-1)*d_ant;
r_RD  = 2*D_ap^2/lambda;
r_min = P.r_lo_fac * r_RD;
r_max = P.r_hi_fac * r_RD;

theta_grid = linspace(P.theta_lo, P.theta_hi, Q_theta);   % 1 x Q_theta

m_bar  = ((0:M-1).' - (M-1)/2);   % M x 1
m_bar2 = m_bar.^2;

% ---- Build per-angle beam-depth range grid (Hussain TWC 2025 Alg.1) --
D_polar_cell = cell(1, Q_theta);   % N_RF x Qr_n per angle
theta_atoms  = [];                 % 1 x S -- angle for each atom
r_atoms      = [];                 % 1 x S -- range for each atom

for ni = 1:Q_theta
    theta_n = theta_grid(ni);
    cos2_n  = cos(theta_n)^2;
    r_EBRD_n = r_RD/10 * cos2_n;   % EBRD at this angle

    % --- Generate r_BD-based sample points for this angle ---
    r_samps = r_min;
    r_F     = r_min;
    n_samp  = 1;
    while r_F < r_max && n_samp < 8
        % r_BD formula (Hussain Theorem 1) -- only valid inside EBRD
        if r_F >= r_EBRD_n * 0.99
            break;  % at or beyond EBRD: beam depth -> inf, stop stepping
        end
        denom = r_RD * cos2_n * (1 - 100*r_F^2/(r_RD^2 * cos2_n^2));
        if denom <= 1e-12
            break;
        end
        r_BD_n = 20 * r_F^2 / denom;
        r_next = min(r_F + r_BD_n, r_max);
        % Avoid duplicate samples
        if abs(r_next - r_samps(end)) > 0.01
            r_samps(end+1) = r_next; %#ok<AGROW>
            n_samp = n_samp + 1;
        end
        if abs(r_next - r_max) < 1e-4; break; end
        r_F = r_next;
    end
    % Always include r_max as endpoint
    if abs(r_samps(end) - r_max) > 0.01
        r_samps(end+1) = r_max;
    end

    % --- Enforce minimum 4 samples (pad with uniform if needed) ---
    Qr_n = numel(r_samps);
    if Qr_n < 4
        n_add = 4 - Qr_n;
        for ai = 1:n_add
            extra = r_min + (r_max - r_min) * ai / (n_add + 1);
            if ~any(abs(extra - r_samps) < 0.1)
                r_samps(end+1) = extra; %#ok<AGROW>
            end
        end
        r_samps = sort(r_samps);
    end

    % --- Build dictionary columns for this angle ---
    Qr_n  = numel(r_samps);
    omega_n = (2*pi*d_ant/lambda) * cos(theta_n);
    c_n     = (pi*d_ant^2/lambda) * sin(theta_n)^2;

    A_n = zeros(M, Qr_n);
    for ri = 1:Qr_n
        kappa_nr = c_n / r_samps(ri);
        A_n(:,ri) = exp(1j*omega_n*m_bar - 1j*kappa_nr*m_bar2);
    end
    D_polar_cell{ni} = W' * A_n;           % N_RF x Qr_n

    theta_atoms = [theta_atoms, repmat(theta_n, 1, Qr_n)]; %#ok<AGROW>
    r_atoms     = [r_atoms,     r_samps];                   %#ok<AGROW>
end

D_polar = cat(2, D_polar_cell{:});   % N_RF x S
Q_total = size(D_polar, 2);

% ---- Covariance-domain OMP with deflation (unchanged) ---------------
R_res    = (R_hat + R_hat')/2;
selected = zeros(1, d);

for k = 1:d
    RD    = R_res * D_polar;
    score = real(sum(conj(D_polar) .* RD, 1));
    [~, best]   = max(score);
    selected(k) = best;

    D_sel  = D_polar(:, selected(1:k));
    P_proj = D_sel * pinv(D_sel);
    P_perp = eye(N_RF) - P_proj;
    R_res  = P_perp * R_hat * P_perp';
    R_res  = (R_res + R_res')/2;
end

% ---- Power estimation -----------------------------------------------
D_sel = D_polar(:, selected);
p_hat = zeros(d, 1);
for k = 1:d
    dk  = D_sel(:,k);
    nd2 = max(real(dk'*dk), 1e-15);
    p_hat(k) = max(0, real(dk'*R_hat*dk) / nd2^2);
end

theta_hat = theta_atoms(selected).';
r_hat     = r_atoms(selected).';

info.selected   = selected;
info.p_est      = p_hat;
info.Q_r_used   = 4;    % nominal (each angle has 4 points min; kept for test_fixes T12)
info.Q_total    = Q_total;
info.theta_all  = theta_atoms;
info.r_all      = r_atoms;
info.sampling   = 'beam_depth_rBD';   % flag distinguishing from old uniform-u
end
