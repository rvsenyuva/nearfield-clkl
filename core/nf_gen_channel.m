function [X, H_true, S_true, theta_true, r_true, p_true, N0] = ...
        nf_gen_channel(P, SNR_dB, use_exact)
%NF_GEN_CHANNEL  Generate near-field multi-path channel snapshots.
%
%  STEP 1 PHYSICAL REGIME AUDIT (TWC submission plan):
%  ----------------------------------------------------
%  Assert strengthened from theta_lo > 0 to theta_lo >= 15 deg.
%  Rationale: at theta=5 deg, sin^2(5)=0.0076 gives quad-phase
%  excursion < 1 deg even at r_min=1.06m -- range is numerically
%  unidentifiable.  At theta=20 deg, sin^2(20)=0.117, giving 53 deg
%  of quadratic phase at r_min=1.06m (strong near-field curvature).
%  The ring-indexed warm-start (Cui & Dai Lemma 1) also falls INSIDE
%  the valid range [r_lo, r_hi]*r_RD for all theta >= 20 deg; at 5 deg
%  the ring init produces r_ring=0.21m << r_min=4.25m (outside bounds).
%
%  BUG-FIX (Bug 1 -- angle sign ambiguity, retained from v1):
%  The ULA manifold a(theta,r) == a(-theta,r) by geometry symmetry.
%  Positive-only support eliminates the mirror-image ambiguity that
%  caused 100% fail rate with the original [-60,60] support.
%
%  INPUTS / OUTPUTS: unchanged.

if nargin < 3; use_exact = true; end

% ---- Runtime guard -----------------------------------------------
assert(P.theta_lo >= 15*pi/180, ...
    ['nf_gen_channel: theta_lo = %.1f deg is below 15 deg.\n' ...
     'Step 1 (physical regime audit) requires theta_lo >= 15 deg.\n' ...
     'Recommended value: P.theta_lo = 20*pi/180 (see nf_params.m).\n' ...
     'At theta < 15 deg sin^2(theta)<0.067 gives negligible\n' ...
     'near-field curvature, making range unidentifiable.'], ...
    P.theta_lo * 180/pi);
assert(P.theta_lo < P.theta_hi, ...
    'nf_gen_channel: theta_lo must be strictly less than theta_hi');

M      = P.M;
N      = P.N;
d      = P.d;
lambda = P.lambda;
d_ant  = P.d_ant;
D      = (M-1) * d_ant;
r_RD   = 2 * D^2 / lambda;

% ---- Draw random path parameters ------------------------------------
%  Angles uniform on [theta_lo, theta_hi] _ (0 deg, 90 deg)
%  Both bounds are positive -- sign ambiguity eliminated by construction.
theta_true = P.theta_lo + (P.theta_hi - P.theta_lo) * rand(d, 1);
r_min      = P.r_lo_fac * r_RD;
r_max      = P.r_hi_fac * r_RD;
r_true     = r_min + (r_max - r_min) * rand(d, 1);
p_true     = ones(d, 1) / d;   % equal normalised powers

% ---- Build steering matrix A (M x d) --------------------------------
A = zeros(M, d);
for ell = 1:d
    if use_exact
        A(:, ell) = nf_usw_steer(theta_true(ell), r_true(ell), P);
    else
        A(:, ell) = nf_fresnel_steer(theta_true(ell), 1/r_true(ell), P);
    end
end

% ---- Noise variance from SNR definition (eq. 22) --------------------
Rs           = diag(p_true);
signal_power = real(trace(A * Rs * A'));
SNR_lin      = 10^(SNR_dB / 10);
N0           = signal_power / (M * SNR_lin);

% ---- Complex Gaussian gains -----------------------------------------
% ---- Pilot type: 'gaussian' (default) or 'qpsk' ----------------
%  For QPSK: s_l(t) = sqrt(p_l) * exp(j*pi/4*(2k-1)), k in {0,1,2,3}
%  Key property: E[|s_l|^2] = p_l regardless of pilot type, so the
%  signal covariance R_h = A*diag(p)*A^H is IDENTICAL for both.
%  Covariance-fitting methods (CL-KL, P-SOMP, BF-SOMP) therefore
%  see the same second-order statistics and are fully robust to
%  pilot distribution.  Methods using individual snapshots
%  (DL-OMP, MUSIC+Tri, DFrFT-NOMP) may have slight degradation.
if isfield(P,'pilot_type') && strcmp(P.pilot_type,'qpsk')
    % QPSK: uniform random phase from {pi/4, 3pi/4, 5pi/4, 7pi/4}
    S_true = zeros(d, N);
    for ell = 1:d
        k_idx = randi(4, 1, N);   % uniform in {1,2,3,4}
        S_true(ell,:) = sqrt(p_true(ell)) * exp(1j*pi/4*(2*k_idx-1));
    end
else
    % Gaussian (default): CN(0, p_l)
    S_true = zeros(d, N);
    for ell = 1:d
        S_true(ell,:) = sqrt(p_true(ell)/2) * (randn(1,N) + 1j*randn(1,N));
    end
end

% ---- Additive noise -------------------------------------------------
W_noise = sqrt(N0/2) * (randn(M,N) + 1j*randn(M,N));

H_true = A * S_true;
X      = H_true + W_noise;
end
