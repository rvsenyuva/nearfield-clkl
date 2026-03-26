function [theta_hat, r_hat, info] = nf_zhang(X, ~, P)
%NF_ZHANG  Dictionary-Learning OMP (DL-OMP) near-field channel estimator.
%
%  Implements Zhang, Zhang & Eldar (IEEE Trans. Commun. 2024), Algorithm 1
%  (multi-path case: Algorithm 2).
%
%  Key idea (_III-A): the distance-parameterised angular dictionary
%     W(r) = [b_N(theta_1,r_1), ..., b_N(theta_N,r_N)]
%  assigns a DIFFERENT r_i to each angle theta_i, making each column unique.
%  This gives much lower inter-column coherence than the 2D polar domain.
%
%  Algorithm (_III-B, Algorithm 2):
%  _________________________________
%  1. Split full array into 2 subarrays of M/2 elements each.
%     Subarray 1: antennas 1..M/2 with reference at centroid_1
%     Subarray 2: antennas M/2+1..M with reference at centroid_2
%     Subarray separation: _ = centroid_2 - centroid_1
%
%  2. Initialise dictionary atoms at distance r_max (Fresnel approx):
%       b_P(theta_i, r_max): standard Fresnel atom, length P/2
%
%  3. for l = 1..L (detect l-th path):
%       for k = 1..K_iter:
%         a. OMP step on subarray 1 residual -> theta__1
%         b. Restrict subarray 2 angles using triangle inequality with
%            known theta__1 and r in [r_min, r_max] -> sub-dictionary
%         c. OMP step on subarray 2 residual -> theta__2
%         d. Distance by law of sines:
%              r_ = _ * sin(theta__2) / sin(theta__2 - theta__1)
%         e. Dictionary update (eq. 23-24):
%              W_j = diag(alpha_j) * W_j^0
%            where alpha_j(p) corrects Fresnel phase from r_max to r_
%       end
%       Update active set, subtract this path from residual
%
%  NOTES ON OUR IMPLEMENTATION:
%  _ X in C^{MxN} is the full (uncompressed) array snapshot matrix.
%    The second argument W (hybrid combiner) is accepted but not used --
%    DL-OMP operates directly on the array (Cui & Dai switching arch.).
%  _ We use K_iter=3 (sufficient per paper's Fig. 9) and L=P.d paths.
%  _ The paper uses a switching architecture selecting P/2 antennas per
%    subarray; we use all M/2 antennas per subarray (equivalent to
%    P = N_RF*T >> M/2, i.e. the full-observation limit, which is valid
%    since X has full dimensions MxN).

M       = P.M;
N       = P.N;
d       = P.d;
lambda  = P.lambda;
d_ant   = P.d_ant;
K_iter  = 3;                   % dictionary update iterations

D_ap  = (M-1)*d_ant;
r_RD  = 2*D_ap^2/lambda;
r_min = P.r_lo_fac * r_RD;
r_max = P.r_hi_fac * r_RD;

% ---- Subarray geometry ---------------------------------------------
M_s    = floor(M/2);
idx1   = 1:M_s;
idx2   = (M-M_s+1):M;
X1     = X(idx1, :);            % M_s x N
X2     = X(idx2, :);            % M_s x N

% Centroid positions along the array axis (0-indexed element centres)
cent1_pos  = mean((idx1 - 1)) * d_ant;
cent2_pos  = mean((idx2 - 1)) * d_ant;
delta      = cent2_pos - cent1_pos;   % subarray separation [m]

% ---- Angle grid and initial Fresnel dictionary at r_max ---------------
Q_th = P.Q_theta;
theta_grid = linspace(P.theta_lo, P.theta_hi, Q_th).';   % Q_th x 1
m_s        = ((0:M_s-1).' - (M_s-1)/2);                   % M_s x 1
m_s2       = m_s.^2;                                       % M_s x 1
omega_g    = (2*pi*d_ant/lambda) * cos(theta_grid).';      % 1 x Q_th

% Initial Fresnel atoms at r_max (will be updated per iteration)
kappa_max = (pi*d_ant^2/lambda) * sin(theta_grid).^2 / r_max;  % Q_th x 1
PHI_max   = m_s * omega_g - m_s2 * kappa_max.';             % M_s x Q_th
W0        = exp(1j * PHI_max);                               % M_s x Q_th (initial dict)

% Per-angle distance estimates (one r per theta, initialised at r_max)
r_active = r_max * ones(d, 1);

% ---- Multi-path OMP loop -------------------------------------------
active_theta = zeros(d, 1);   % estimated angles (rad)
active_r     = zeros(d, 1);   % estimated ranges (m)

% Build subarray covariances
R1_full = (1/N) * (X1 * X1');   % M_s x M_s
R2_full = (1/N) * (X2 * X2');

z1 = R1_full;   % residual covariance for subarray 1
z2 = R2_full;   % residual covariance for subarray 2

for l = 1:d

    % Current angle-specific dictionaries (updated each iteration)
    W1 = W0;   W2 = W0;

    theta1_hat = theta_grid(1);
    theta2_hat = theta_grid(end);
    r_hat_l    = 0.5*(r_min + r_max);

    for k = 1:K_iter
        % ---- Step 1: angle estimate from subarray 1 -----------------
        sc1 = real(sum(conj(W1) .* (z1 * W1), 1)).';  % Q_th x 1
        [~, idx1_best] = max(sc1);
        theta1_hat = theta_grid(idx1_best);

        % ---- Step 2: constrained angle search in subarray 2 ---------
        % From theta__1 and r in [r_min, r_max], the angle seen by subarray 2
        % is bounded (eq. 20 in the paper).  Approximate via angle range
        % [theta_lo, theta_hi] (same grid -- triangulation disambiguates).
        % For simplicity we search the full grid and rely on the
        % triangulation step to produce a self-consistent solution.
        sc2 = real(sum(conj(W2) .* (z2 * W2), 1)).';  % Q_th x 1
        [~, idx2_best] = max(sc2);
        theta2_hat = theta_grid(idx2_best);

        % Avoid degenerate case where both subarrays report same angle
        if abs(theta2_hat - theta1_hat) < 1e-4
            % Perturb search window for subarray 2 to second-best peak
            sc2_tmp = sc2;
            nw = max(1, round(Q_th/30));
            sc2_tmp(max(1,idx2_best-nw):min(Q_th,idx2_best+nw)) = -Inf;
            [~, idx2_best] = max(sc2_tmp);
            theta2_hat = theta_grid(idx2_best);
        end

        % ---- Step 3: distance by law of sines -----------------------
        denom_r = sin(theta2_hat) - sin(theta1_hat);
        if abs(denom_r) < 1e-6
            r_hat_l = 0.5*(r_min + r_max);
        else
            r_hat_l = delta * sin(theta2_hat) / ...
                      max(1e-6, abs(theta2_hat - theta1_hat)) * ...
                      sign(sin(theta2_hat - theta1_hat));
            % Use law of sines formula from eq. (18):
            % r1/sin(theta2) = _/sin(theta2-theta1)
            dsin = sin(theta2_hat - theta1_hat);
            if abs(dsin) > 1e-6
                r_hat_l = delta * sin(theta2_hat) / dsin;
            end
        end
        r_hat_l = min(r_max, max(r_min, abs(r_hat_l)));

        % ---- Step 4: dictionary update (eq. 23-24) ------------------
        %  alpha(p) = exp(-j*(p-1)^2*(_ d^2 sin^2theta_/lambda) * (1/r_ - 1/r_max))
        %  This corrects the phase of each atom from r_max to r_hat_l
        u_hat  = 1/r_hat_l;
        u_max  = 1/r_max;

        % Update subarray 1 dictionary (all angles, estimated theta__1)
        correction1 = exp(-1j * m_s2 * ((pi*d_ant^2/lambda)*sin(theta1_hat)^2*(u_hat-u_max)));
        W1(:, idx1_best) = correction1 .* W0(:, idx1_best);

        % Update subarray 2 dictionary (all angles, estimated theta__2)
        correction2 = exp(-1j * m_s2 * ((pi*d_ant^2/lambda)*sin(theta2_hat)^2*(u_hat-u_max)));
        W2(:, idx2_best) = correction2 .* W0(:, idx2_best);
    end

    % Averaged angle and range for this path
    active_theta(l) = (theta1_hat + theta2_hat) / 2;
    active_r(l)     = r_hat_l;

    % Orthogonal projection deflation on both subarrays
    a1_l = W1(:, idx1_best);
    a2_l = W2(:, idx2_best);

    P1     = a1_l * pinv(a1_l);
    Perp1  = eye(M_s) - P1;
    z1     = Perp1 * R1_full * Perp1';
    z1     = (z1 + z1')/2;

    P2     = a2_l * pinv(a2_l);
    Perp2  = eye(M_s) - P2;
    z2     = Perp2 * R2_full * Perp2';
    z2     = (z2 + z2')/2;
end

theta_hat = active_theta;
r_hat     = active_r;

info.theta1_per_path = active_theta;
info.r_per_path      = active_r;
info.K_iter          = K_iter;
info.delta_m         = delta;
end
