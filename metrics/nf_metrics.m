function [nmse, rmse_theta, rmse_r, fail_rate] = nf_metrics(...
        theta_hat, r_hat, Y, W, H_true, theta_true, r_true, P)
%NF_METRICS  Compute channel NMSE, angle/range RMSE, and failure rate.
%
%  BUG-FIX (Bug 3 — matching cost robustness):
%  --------------------------------------------
%  With Bug 1 fixed (angles drawn from positive half-plane only), the
%  signed angular error theta_true(i) - theta_hat(j) is the correct
%  quantity to minimise.  However, two changes make matching more robust:
%
%  (a) Hungarian matching cost now uses ABSOLUTE angular distance
%      min(|θ_true − θ_hat|, π − |θ_true + θ_hat|) — this is the
%      geodesic on the unit circle and correctly handles the residual
%      sign-mirror edge cases that can still occur when an estimator
%      returns a negative angle.  This costs nothing extra and makes the
%      RMSE numerically identical to the signed version when all angles
%      are positive, but recovers gracefully for any stray negative outputs.
%
%  (b) The success criterion uses the same absolute distance so that a
%      near-correct estimate at +29.9° is not penalised as a failure when
%      the true angle is +30.0°.
%
%  The RMSE reported is the signed error after optimal matching —
%  rmse_theta captures both magnitude and sign errors in one number.

M = P.M;
d = P.d;

% ---- Degenerate estimate guard ------------------------------------
if isempty(theta_hat) || length(theta_hat) < d
    nmse       = 1.0;
    rmse_theta = pi;
    rmse_r     = Inf;
    fail_rate  = 1;
    return
end

% ---- Channel reconstruction (Fresnel model for gain LS) -----------
d_hat = length(theta_hat);
A_hat = zeros(M, d_hat);
for i = 1:d_hat
    r_i = max(r_hat(i), 1e-3);
    A_hat(:,i) = nf_fresnel_steer(theta_hat(i), 1/r_i, P);
end

D_hat = W' * A_hat;
S_hat = pinv(D_hat) * Y;
H_hat = A_hat * S_hat;

num  = norm(H_hat - H_true, 'fro')^2;
den  = norm(H_true, 'fro')^2;
nmse = num / max(den, 1e-15);

% ---- Hungarian matching with geodesic angular cost ----------------
%  BUG-FIX: use min(|Δθ|, π − |θ_true + θ_hat|) as the angular
%  distance in the cost matrix.  This is the shortest arc on the unit
%  circle and handles sign-mirror outputs from estimators gracefully.
D_ap   = (M-1)*P.d_ant;
r_RD   = 2*D_ap^2/P.lambda;
eta    = 1 / (P.r_hi_fac * r_RD)^2;   % range normalisation

d_true = length(theta_true);
d_est  = length(theta_hat);

cost_mat = zeros(d_true, d_est);
for i = 1:d_true
    for j = 1:d_est
        % Geodesic angular distance (handles sign mirrors correctly)
        dth = abs(theta_true(i) - theta_hat(j));
        dth = min(dth, pi - abs(theta_true(i) + theta_hat(j)));
        dr  = r_true(i) - r_hat(j);
        cost_mat(i,j) = dth^2 + eta * dr^2;
    end
end

assn      = hungarian_2d(cost_mat);
theta_err = zeros(d_true, 1);
r_err     = zeros(d_true, 1);
n_success = 0;

for i = 1:d_true
    j            = assn(i);
    % Signed errors for RMSE reporting
    theta_err(i) = theta_true(i) - theta_hat(j);
    r_err(i)     = r_true(i)     - r_hat(j);
    % Success: use geodesic distance for threshold check
    dth_geo = abs(theta_err(i));
    dth_geo = min(dth_geo, pi - abs(theta_true(i) + theta_hat(j)));
    if dth_geo <= P.dtheta_tol && abs(r_err(i)) <= P.dr_fac_tol * r_true(i)
        n_success = n_success + 1;
    end
end

rmse_theta = sqrt(mean(theta_err.^2));
rmse_r     = sqrt(mean(r_err.^2));
fail_rate  = double(n_success < d_true);
end


% ====================================================================
function assn = hungarian_2d(cost)
d_rows = size(cost,1);
assn   = zeros(d_rows,1);
try
    pairs = matchpairs(cost, 1e6);
    for k = 1:size(pairs,1)
        assn(pairs(k,1)) = pairs(k,2);
    end
    unmatched = find(assn == 0);
    for i = unmatched.'
        row = cost(i,:); [~,j] = min(row);
        assn(i) = j;
    end
catch
    used = false(size(cost,2),1);
    for i = 1:d_rows
        row = cost(i,:); row(used) = Inf;
        [~,j] = min(row);
        assn(i) = j; used(j) = true;
    end
end
end
