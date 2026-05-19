function [theta_hat, r_hat, p_hat, u_hat, info] = nf_fckl(R_hat, W, P)
%NF_FCKL  Fixed-Curvature KL estimator (ablation-style 7th baseline).
%
%  Thin wrapper around nf_clkl with do_post_loop_scan = false.
%  This ablates the post-loop global matched-filter scan (Phase 2),
%  so curvature estimates are frozen at the warm-start values found by
%  the multi-start OMP initialisation.  The power-only main loop
%  (Phase 1) runs normally with the best warm-start support.
%
%  Used in Fig. 2 and related sweeps as a compressed-domain peer that
%  shares the same KL objective but lacks Phase 2 curvature refinement.
%
%  INPUTS / OUTPUTS: identical to nf_clkl.m (see that file for docs).
%
%  See also: nf_clkl

[theta_hat, r_hat, p_hat, u_hat, info] = nf_clkl(R_hat, W, P, false);
end
