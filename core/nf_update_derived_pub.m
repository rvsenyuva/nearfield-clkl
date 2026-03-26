function P = nf_update_derived_pub(P)
%NF_UPDATE_DERIVED_PUB  Public wrapper — recompute aperture D and Rayleigh
%  distance r_RD after changing P.M.  Also resets u bounds.
%
%  P = nf_update_derived_pub(P)
%
%  Call this any time P.M is changed before running a simulation.

P.D    = (P.M - 1) * P.d_ant;
P.r_RD = 2 * P.D^2 / P.lambda;

% Recompute u bounds from current r_lo/hi factors + margin
r_min_sim = P.r_lo_fac * P.r_RD / P.u_margin;
r_max_sim = P.r_hi_fac * P.r_RD * P.u_margin;
P.u_min   = 1 / r_max_sim;
P.u_max   = 1 / r_min_sim;
end
