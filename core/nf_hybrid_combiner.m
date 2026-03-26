function [W, Y, R_hat] = nf_hybrid_combiner(X, P)
%NF_HYBRID_COMBINER  Build random-phase combiner, compress, whiten, form sample covariance.
%
%  [W, Y, R_hat] = nf_hybrid_combiner(X, P)
%
%  Implements the hybrid observation model (eq. 2):
%      y(n) = W^H * x(n)  in C^{N_RF}
%
%  Combiner design (eq. 16):
%      [W]_{m,k} = (1/sqrt(M)) * exp(j*phi_{m,k}),   phi ~ U[0, 2*pi)
%
%  Whitening (to make effective noise spatially white in compressed domain):
%      y_tilde(n) = (W^H*W)^{-1/2} * y(n)
%
%  Sample covariance (eq. 3):
%      R_hat = (1/N) * Y_tilde * Y_tilde^H
%
%  INPUTS
%    X  : M x N snapshot matrix
%    P  : parameter struct (uses P.M, P.N, P.N_RF)
%
%  OUTPUTS
%    W      : M x N_RF  constant-modulus combining matrix
%    Y      : N_RF x N  whitened compressed snapshots
%    R_hat  : N_RF x N_RF  sample covariance of whitened snapshots

M    = P.M;
N    = P.N;
N_RF = P.N_RF;

% ---- Constant-modulus (random-phase) combiner -----------------------
phi  = 2*pi * rand(M, N_RF);
W    = (1/sqrt(M)) * exp(1j * phi);          % M x N_RF

% ---- Compressed snapshots -------------------------------------------
Y_raw = W' * X;                              % N_RF x N

% ---- Whitening  (W^H*W may differ from I for random W) --------------
WtW  = W' * W;                               % N_RF x N_RF
% Symmetric positive-definite square root inverse via eigendecomposition
[V, D_eig] = eig((WtW + WtW')/2);
d_eig      = real(diag(D_eig));
d_eig      = max(d_eig, 1e-12);             % numerical floor
WtW_inv_sqrt = V * diag(1./sqrt(d_eig)) * V';

Y    = WtW_inv_sqrt * Y_raw;                 % N_RF x N  (whitened)

% ---- Sample covariance ----------------------------------------------
R_hat = (1/N) * (Y * Y');                   % N_RF x N_RF
R_hat = (R_hat + R_hat') / 2;               % enforce Hermitian symmetry
end
