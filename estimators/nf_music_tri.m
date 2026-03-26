function [theta_hat, r_hat, info] = nf_music_tri(X, W_or_P, P_arg)
%NF_MUSIC_TRI  Subarray MUSIC + geometric triangulation.
%
%  Calling conventions (both supported):
%    [th,r,info] = nf_music_tri(X, P)        -- original 2-arg form
%    [th,r,info] = nf_music_tri(X, W, P)     -- 3-arg form (W ignored;
%                                                MUSIC+Tri operates on
%                                                full-array X subarrays
%                                                directly and does not use
%                                                the hybrid combiner W)
%
%  CHANGES in this version:
%  -----------------------------------------------------------------------
%  1. ROBUST LS TRIANGULATION: The Q=3 LS system (Haghshenas eq. 17-18)
%     can be ill-conditioned when two sources are nearly co-linear with
%     a subarray baseline.  We now check rcond(A_k'*A_k) and fall back
%     to the pairwise law-of-sines (Q=2 equivalent) for any source whose
%     LS system is ill-conditioned (rcond < 1e-8).
%
%  2. MATCHED ANGLE SIGN: all subarray angle estimates are clipped to
%     [P.theta_lo, P.theta_hi] so negative-angle artefacts cannot
%     propagate into the triangulation.
%
%  3. SCAN RANGE RESTRICTION (v5 bug-fix retained): scan only on
%     [P.theta_lo, P.theta_hi] to prevent mirror-image peaks.
%
%  4. Q = P.n_subarrays (default 3).  For ill-conditioned Q=3 cases the
%     fallback produces a 2-baseline estimate instead of NaN/Inf.

% ---- Argument dispatch: support (X,P) and (X,W,P) -------------------
if nargin == 2
    P = W_or_P;   % old-style call: nf_music_tri(X, P)
else
    P = P_arg;    % new-style call: nf_music_tri(X, W, P)  [W not used]
end

M     = P.M;
N     = P.N;
d     = P.d;
d_ant = P.d_ant;

Q     = P.n_subarrays;
M_s   = floor(M / Q);

idx_sub = cell(Q, 1);
R_sub   = cell(Q, 1);
pos_sub = zeros(Q, 1);

for q = 1:Q
    idx_sub{q} = (q-1)*M_s + (1:M_s);
    X_q        = X(idx_sub{q}, :);
    R_sub{q}   = (1/N) * (X_q * X_q');
    pos_sub(q) = mean(idx_sub{q} - 1) * d_ant;
end

% ---- MUSIC angle estimation per subarray ----------------------------
theta_sub = zeros(d, Q);
for q = 1:Q
    theta_sub(:, q) = music_spectrum_fast(R_sub{q}, d, M_s, d_ant, ...
                                          P.lambda, P.theta_lo, P.theta_hi);
end

% Sort ascending in each subarray (Algorithm 1 step 6 of Haghshenas)
for q = 1:Q
    theta_sub(:, q) = sort(theta_sub(:, q));
end

% Cross-subarray Hungarian matching
cost12 = abs(bsxfun(@minus, theta_sub(:,1), theta_sub(:,2).'));
assn12 = hungarian_assign(cost12);
theta_sub(:, 2) = theta_sub(assn12, 2);

if Q >= 3
    cost13 = abs(bsxfun(@minus, theta_sub(:,1), theta_sub(:,3).'));
    assn13 = hungarian_assign(cost13);
    theta_sub(:, 3) = theta_sub(assn13, 3);
end

% Clip all angles to valid range (prevents negative-angle artefacts)
theta_sub = max(P.theta_lo, min(P.theta_hi, theta_sub));

% ---- Triangulation --------------------------------------------------
D_ap  = (M-1)*d_ant;
r_RD  = 2*D_ap^2/P.lambda;
r_fallback = 0.5*(P.r_lo_fac + P.r_hi_fac) * r_RD;

theta_hat = zeros(d, 1);
r_hat     = zeros(d, 1);

for k = 1:d
    if Q == 2
        % Direct cotangent formula
        denom = cot(theta_sub(k,1)) - cot(theta_sub(k,2));
        if abs(denom) < 1e-6
            r_hat(k) = r_fallback;
        else
            r_hat(k) = abs((pos_sub(2) - pos_sub(1)) / denom);
        end
        theta_hat(k) = mean(theta_sub(k,:));
    else
        % Full LS triangulation (Haghshenas eq. 17-18) with rcond fallback
        A_k = zeros(2*(Q-1), Q);
        b_k = zeros(2*(Q-1), 1);
        p1  = [0; pos_sub(1)];

        for qi = 2:Q
            rows = (qi-2)*2 + (1:2);
            pq   = [0; pos_sub(qi)];
            dk1  = [sin(theta_sub(k,1)); cos(theta_sub(k,1))];
            dkq  = [sin(theta_sub(k,qi)); cos(theta_sub(k,qi))];
            b_k(rows)     = p1 - pq;
            A_k(rows, 1)  = -dk1;
            A_k(rows, qi) =  dkq;
        end

        AtA = A_k.' * A_k;
        % Use rcond to detect ill-conditioned systems
        if rcond(AtA) > 1e-10
            t_hat_k = AtA \ (A_k.' * b_k);
            uk_estimates = zeros(2, Q);
            for qi = 1:Q
                pq  = [0; pos_sub(qi)];
                dkq = [sin(theta_sub(k,qi)); cos(theta_sub(k,qi))];
                uk_estimates(:, qi) = pq + t_hat_k(qi) * dkq;
            end
            uk_mean      = mean(uk_estimates, 2);
            r_hat(k)     = norm(uk_mean);
            theta_hat(k) = atan2(uk_mean(1), uk_mean(2));
            theta_hat(k) = abs(theta_hat(k));
        else
            % Fallback: use first two subarrays only (law of sines)
            denom = cot(theta_sub(k,1)) - cot(theta_sub(k,2));
            if abs(denom) < 1e-6
                r_hat(k) = r_fallback;
            else
                r_hat(k) = abs((pos_sub(2) - pos_sub(1)) / denom);
            end
            theta_hat(k) = mean(theta_sub(k,:));
        end
    end
end

% Clip to physically valid range
r_hat = max(P.r_lo_fac*r_RD/2, min(P.r_hi_fac*r_RD*2, r_hat));
% Clip theta to valid range
theta_hat = max(P.theta_lo, min(P.theta_hi, theta_hat));

info.theta_sub = theta_sub;
info.Q_used    = Q;
end


% ====================================================================
function peaks = music_spectrum_fast(R, d, M_s, d_ant, lambda, theta_lo, theta_hi)
N_scan     = 512;
theta_scan = linspace(theta_lo, theta_hi, N_scan);

[V, D_eig] = eig((R+R')/2);
ev         = real(diag(D_eig));
[~, ord]   = sort(ev, 'ascend');
n_noise    = max(1, M_s - d);
E_n        = V(:, ord(1:n_noise));

m_bar_s = ((0:M_s-1).' - (M_s-1)/2);
omega_s = (2*pi*d_ant/lambda) * cos(theta_scan);
A_scan  = exp(1j * m_bar_s * omega_s);

EtA   = E_n' * A_scan;
pspec = 1 ./ max(sum(abs(EtA).^2, 1), 1e-15);

peaks = zeros(d, 1);
pspec_local = pspec;
nw = max(1, round(N_scan/50));
for k = 1:d
    [~, idx]  = max(pspec_local);
    peaks(k)  = theta_scan(idx);
    pspec_local(max(1,idx-nw):min(N_scan,idx+nw)) = 0;
end
end


function assn = hungarian_assign(cost)
d    = size(cost, 1);
assn = zeros(d, 1);
used = false(d, 1);
for i = 1:d
    row = cost(i,:); row(used) = Inf;
    [~, j] = min(row);
    assn(i) = j; used(j) = true;
end
end
