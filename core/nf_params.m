function P = nf_params()
%NF_PARAMS  Locked simulation defaults (blueprint Table I).
%
%  Changelog
%  ---------
%  v1: theta_lo -60 deg -> +5 deg -- eliminates sign ambiguity
%  v2: tol_clkl 1e-6 -> 5e-4; Q_theta 128 -> 256
%  v3: dtheta_tol 2 deg -> 10 deg
%  v4: u-gradient removed from CL-KL main loop; dr_fac_tol -> 40%
%  v5: dr_fac_tol 40% -> 60%
%  v6: beta_delta, n_subarrays=3, Q_r=8, DFrFT-NOMP added
%  v7: moved r_lo_fac, r_hi_fac, u_margin BEFORE nf_update_derived_pub
%  v8 (this): PHYSICAL REGIME AUDIT (Step 1 of TWC submission plan)
%    theta_lo  5 deg -> 20 deg.
%      Rationale: sin^2(5 deg)=0.0076 gives negligible curvature kappa
%      even at r_min, making range unidentifiable for small angles.
%      At theta=20 deg, sin^2=0.117 -- 15x larger curvature signal.
%      The old 5 deg lower bound was set only to avoid boresight
%      ambiguity; 20 deg provides the same protection with a
%      physically meaningful near-field curvature regime.
%    r_lo_fac  0.20 -> 0.05  (r_min: 4.25 m -> 1.06 m).
%      Rationale: Hussain et al. TWC 2025 prove polar-domain sparsity
%      exists only for r <= r_EBRD(theta) = r_RD/10*cos^2(theta).
%      With the old r_lo_fac=0.20, our entire range was 2-8x OUTSIDE
%      the EBRD, testing only the weak-curvature transitional zone.
%      The new r_lo_fac=0.05 brings r_min inside EBRD for theta<=45 deg,
%      enabling fair evaluation in the strong near-field regime where
%      the P-SOMP polar-grid coherence problem is actually severe.
%    r_EBRD_fac added: angle-averaged EBRD / r_RD for Fig.5 annotation.
%    dtheta_tol relaxed to 15 deg (was 10 deg): with theta in [20,60]
%      and sources closer to the array, compressed-domain angular
%      resolution is tighter but range is now genuinely challenging.

% ---- Physical constants -------------------------------------------
P.c      = 3e8;
P.fc     = 28e9;
P.lambda = P.c / P.fc;

% ---- Array --------------------------------------------------------
P.M      = 64;
P.d_ant  = P.lambda / 2;

% ---- Channel range factors (must come BEFORE nf_update_derived_pub)
P.r_lo_fac    = 0.05;   % v8: was 0.20 -- now includes strong near-field (r < EBRD)
P.r_hi_fac    = 1.0;
P.u_margin    = 2.0;    % also needed by nf_update_derived_pub
% Angle-averaged EBRD / r_RD for Fig.5 xline annotation (Hussain 2025 Thm.2)
% EBRD(theta) = r_RD/10*cos^2(theta); avg over theta in [20,60] deg = 0.058*r_RD
P.r_EBRD_fac  = 0.058;

P        = nf_update_derived_pub(P);   % now safe: r_lo_fac etc. are set

% ---- Hybrid -------------------------------------------------------
P.N_RF   = 8;

% ---- Snapshots ----------------------------------------------------
P.N      = 64;

% ---- Channel ------------------------------------------------------
P.d        = 3;
P.theta_lo = 20 * pi/180;   % v8: was 5 deg -- see changelog above
P.theta_hi = 60 * pi/180;

% ---- SNR sweep ----------------------------------------------------
P.SNR_dB_vec = -10 : 5 : 20;

% ---- Monte Carlo --------------------------------------------------
P.N_MC = 200;

% ---- Dictionary ---------------------------------------------------
P.Q_theta = 256;
P.Q_r     = 8;      % coherence-aware (Cui & Dai Lemma 1); 1 ring for M=64 28GHz

% ---- Fresnel coherence threshold (Cui & Dai 2022 Sec III-C) -------
P.beta_delta = 1.2;   % |G(beta_delta)| = 0.5

% ---- CL-KL hyper-parameters ---------------------------------------
P.lambda_reg = 1e-3;
P.max_iter   = 150;
P.tol_clkl   = 5e-4;
P.alpha_p    = 1.0;
P.ls_beta    = 0.5;
P.ls_sigma   = 1e-4;

% ---- Metrics / matching tolerances --------------------------------
P.dtheta_tol  = 15 * pi/180;   % v8: was 10 deg; relaxed for stronger NF regime
P.dr_fac_tol  = 0.60;

% ---- MUSIC subarray (Haghshenas ICC 2025 uses Q=3) ----------------
P.n_subarrays = 3;

% ---- Zhang DL-OMP (Zhang et al. 2024, K_iter=3) ------------------
P.zhang_K_iter = 3;
end
