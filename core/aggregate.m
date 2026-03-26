function [nmse_pt, nstd_pt, rth_pt, rr_pt, fail_pt, ...
          clkl_iters_mean, clkl_conv_pct, clkl_N0_mean] = ...
    aggregate(nmse_mc, rmse_th_mc, rmse_r_mc, fail_mc, n_meth, ...
              clkl_iters_mc, clkl_conv_mc, clkl_N0_mc)
%AGGREGATE  Compute mean/std across Monte Carlo trials for all metrics.
%
%  [nmse_pt, nstd_pt, rth_pt, rr_pt, fail_pt, ...
%   clkl_iters_mean, clkl_conv_pct, clkl_N0_mean] = ...
%       aggregate(nmse_mc, rmse_th_mc, rmse_r_mc, fail_mc, n_meth, ...
%                 clkl_iters_mc, clkl_conv_mc, clkl_N0_mc)
%
%  Required inputs (all n_meth x N_MC matrices):
%    nmse_mc      — raw linear NMSE per trial (converted to dB here)
%    rmse_th_mc   — angle RMSE in degrees per trial
%    rmse_r_mc    — range RMSE in metres per trial
%    fail_mc      — failure indicator (0 or 100) per trial
%    n_meth       — number of methods
%
%  Optional CL-KL diagnostic inputs (1 x N_MC row vectors, pass [] to skip):
%    clkl_iters_mc — n_iter per trial for CL-KL (method 1)
%    clkl_conv_mc  — converged flag (0/1) per trial for CL-KL
%    clkl_N0_mc    — final N0_hat per trial for CL-KL
%
%  New diagnostic outputs (scalars, NaN if not provided):
%    clkl_iters_mean — mean CL-KL iterations across trials
%    clkl_conv_pct   — % of trials where CL-KL converged (hit tol, not max_iter)
%    clkl_N0_mean    — mean final N0 estimate (collapse indicator: <1e-6 = bad)

if nargin < 6; clkl_iters_mc = []; end
if nargin < 7; clkl_conv_mc  = []; end
if nargin < 8; clkl_N0_mc    = []; end

nmse_pt = nan(n_meth,1);
nstd_pt = nan(n_meth,1);
rth_pt  = nan(n_meth,1);
rr_pt   = nan(n_meth,1);
fail_pt = nan(n_meth,1);

for mi = 1:n_meth
    v           = 10*log10(nmse_mc(mi,:));
    nmse_pt(mi) = mean(v, 'omitnan');
    nstd_pt(mi) = std(v,  'omitnan');
    rth_pt(mi)  = mean(rmse_th_mc(mi,:), 'omitnan');
    rr_pt(mi)   = mean(rmse_r_mc(mi,:),  'omitnan');
    fail_pt(mi) = mean(fail_mc(mi,:),    'omitnan');
end

% CL-KL diagnostics
if ~isempty(clkl_iters_mc)
    clkl_iters_mean = mean(clkl_iters_mc, 'omitnan');
else
    clkl_iters_mean = NaN;
end

if ~isempty(clkl_conv_mc)
    clkl_conv_pct = 100 * mean(clkl_conv_mc, 'omitnan');
else
    clkl_conv_pct = NaN;
end

if ~isempty(clkl_N0_mc)
    clkl_N0_mean = mean(clkl_N0_mc, 'omitnan');
else
    clkl_N0_mean = NaN;
end
end
