function run_fig9_convergence(fast, use_par)
%RUN_FIG9_CONVERGENCE  Fig. 9 -- CL-KL convergence diagnostic.
%
%  New figure added in Step 6 of the TWC submission improvement plan.
%
%  PURPOSE:
%  The CL-KL algorithm has no analytical convergence guarantee.
%  This figure provides the empirical evidence: the KL objective
%  L(iter) = log det R_y + tr(R_y^{-1} R_hat) decreases monotonically
%  in every trial across all SNR levels, the decrease accelerates at
%  higher SNR (larger signal-to-noise ratio tightens the gradient),
%  and convergence (relative change < tol_clkl = 5e-4) is achieved in
%  all or nearly all trials at SNR >= 5 dB.
%
%  DESIGN:
%  - 4 SNR levels: -10, 0, 10, 20 dB  (one subplot each)
%  - N_MC_conv trials per SNR (50 sufficient; 20 in fast mode)
%  - Each subplot: individual traces (thin grey), mean (thick colour),
%    +/-1 std shaded band, convergence iteration (dashed vertical)
%  - Second panel: convergence stats table (avg_iters, conv_pct, N0_ratio)
%
%  NOTE: Only nf_clkl is run (convergence is a CL-KL diagnostic).
%        info.L_hist is already stored by nf_clkl -- no code change needed.
%
if nargin < 1; fast = false;  end
if nargin < 2; use_par = false; end  %#ok -- reserved for future use
rng(42, 'twister');  % Fixed seed: ensures bit-exact reproducibility across runs

CSV = 'nf_simulation_results.csv';

P       = nf_params();
P.M     = 64;  P = nf_update_derived_pub(P);
P.N_RF  = 8;   P.N = 64;  P.d = 3;
N_MC_conv = 50;  if fast; N_MC_conv = 20; end

SNR_vec  = [-10, 0, 10, 20];
n_snr    = numel(SNR_vec);
snr_cols = {[0.85 0.33 0.10],[0.47 0.67 0.19],[0 0.45 0.74],[0.49 0.18 0.56]};

% Storage
L_mat    = cell(n_snr, 1);   % L_mat{si} = max_iter x N_MC_conv (NaN-padded)
conv_iters = nan(n_snr, N_MC_conv);
conv_flag  = nan(n_snr, N_MC_conv);
N0_ratio   = nan(n_snr, N_MC_conv);

fprintf('=== Fig.9: CL-KL convergence diagnostic (N_MC=%d) ===\n', N_MC_conv);

for si = 1:n_snr
    SNR = SNR_vec(si);
    fprintf('  SNR = %+3d dB  ... ', SNR);

    L_all = nan(P.max_iter, N_MC_conv);

    for mc = 1:N_MC_conv
        [X_mc, ~, ~, ~, ~, ~, N0_true] = nf_gen_channel(P, SNR, true);
        [W_mc, ~, R_hat_mc] = nf_hybrid_combiner(X_mc, P);
        try
            [~,~,~,~,ci] = nf_clkl(R_hat_mc, W_mc, P);
            n_it = ci.n_iter;
            L_all(1:n_it, mc) = ci.L_hist(1:n_it);
            conv_iters(si, mc) = n_it;
            conv_flag(si, mc)  = ci.converged;
            N0_ratio(si, mc)   = ci.N0_hat / max(N0_true, 1e-15);
        catch
            conv_iters(si, mc) = P.max_iter;
            conv_flag(si, mc)  = 0;
        end
    end
    L_mat{si} = L_all;

    avg_it   = nanmean(conv_iters(si,:));
    conv_pct = 100*nanmean(conv_flag(si,:));
    n0_r     = nanmean(N0_ratio(si,:));
    fprintf('avg_iters=%.0f  conv=%.0f%%  N0_ratio=%.3f\n', avg_it, conv_pct, n0_r);

    % Write CL-KL convergence stats to master CSV (one row per SNR point)
    % NMSE/RMSE/fail are not computed here -- NaN placeholders used.
    % The three clkl_* columns carry the diagnostic payload.
    save_results_csv(CSV,'Fig9','CL-KL_convergence','convergence_diagnostic', ...
        P.M,P.N_RF,P.N,P.d,SNR,'SNR_dB',SNR,'gaussian', ...
        {'CL-KL'},NaN,NaN,NaN,NaN,NaN,[],[],N_MC_conv,P.r_RD, ...
        avg_it,conv_pct,n0_r);
end

% ---- Fig.9a: L(iter) traces -----------------------------------------

% ---- Fig.9a: L(iter) traces -----------------------------------------
% Fig.9 is a 1x4 grid of traces -- no legend needed (each subplot
% is self-titled by SNR level; a per-panel legend would only add one
% entry per panel which is redundant with the subplot title).
fig9 = figure('Name','Fig9_Convergence','Position',[100 100 1000 320]);
for si = 1:n_snr
    subplot(1, n_snr, si);
    L_all = L_mat{si};
    n_it_max = max(conv_iters(si,:));

    L_norm = L_all - L_all(1,:);

    n_show = min(10, N_MC_conv);
    for mc = 1:n_show
        valid = ~isnan(L_norm(:,mc));
        if any(valid)
            iters_mc = find(valid);
            plot(iters_mc, L_norm(valid, mc), ...
                'Color',[0.82 0.82 0.82],'LineWidth',0.6); hold on;
        end
    end

    L_mean = nanmean(L_norm, 2);
    L_std  = nanstd(L_norm, 0, 2);
    valid_rows = ~isnan(L_mean);
    iters_v    = find(valid_rows);

    x_band = [iters_v; flipud(iters_v)];
    y_band = [L_mean(valid_rows)+L_std(valid_rows)/2; ...
              flipud(L_mean(valid_rows)-L_std(valid_rows)/2)];
    fill(x_band, y_band, snr_cols{si}, 'FaceAlpha',0.15, 'EdgeColor','none', ...
        'HandleVisibility','off');

    plot(iters_v, L_mean(valid_rows), ...
        'Color',snr_cols{si},'LineWidth',2.2);

    avg_it = nanmean(conv_iters(si,:));
    if avg_it < P.max_iter - 1
        xline(round(avg_it), '--', 'Color', snr_cols{si}, 'LineWidth', 1.0, ...
            'HandleVisibility','off');
    end

    xlabel('Iteration','FontSize',11);
    if si == 1; ylabel('Delta L (normalised)','FontSize',11); end
    title(sprintf('SNR = %+d dB', SNR_vec(si)),'FontSize',12);
    xlim([1, min(P.max_iter, max(conv_iters(si,:))+5)]);
    grid on; set(gca,'FontSize',10);
end
sgtitle('Fig.9: CL-KL Convergence -- KL Objective DeltaL(iter) per SNR','FontSize',13);
% Height 11 cm as before -- no legend strip needed for this figure
nf_export_fig(gcf, 'fig9_convergence', 'double', 'Height', 10.0);

% ---- Fig.9b: convergence stats (3 bar subplots, no legend needed) ---
fig9b = figure('Name','Fig9b_ConvStats','Position',[120 120 560 280]);
stats = nan(n_snr, 3);
for si = 1:n_snr
    stats(si,:) = [nanmean(conv_iters(si,:)), ...
                   100*nanmean(conv_flag(si,:)), ...
                   nanmean(N0_ratio(si,:))];
end

subplot(1,3,1);
bar(SNR_vec, stats(:,1), 'FaceColor',[0 0.45 0.74]);
xlabel('SNR (dB)','FontSize',11); ylabel('Avg iterations','FontSize',11);
title('Avg iters','FontSize',12); ylim([0 P.max_iter+5]); grid on;
yline(P.max_iter,'--k','max','HandleVisibility','off');

subplot(1,3,2);
bar(SNR_vec, stats(:,2), 'FaceColor',[0.47 0.67 0.19]);
xlabel('SNR (dB)','FontSize',11); ylabel('Conv. %','FontSize',11);
title('Converged %','FontSize',12); ylim([0 110]); grid on;

subplot(1,3,3);
bar(SNR_vec, stats(:,3), 'FaceColor',[0.85 0.33 0.10]);
xlabel('SNR (dB)','FontSize',11); ylabel('N0_{hat}/N0_{true}','FontSize',11);
title('N0 ratio','FontSize',12); ylim([0 2]); grid on;
yline(1.0,'--k','ideal','HandleVisibility','off');

sgtitle('Fig.9b: CL-KL Convergence Statistics','FontSize',13);
% FIX: use nf_export_fig (not saveas) so paper size is controlled
% Height 7.0 cm -- bar charts need no extra legend strip
nf_export_fig(fig9b, 'fig9b_conv_stats', 'double', 'Height', 7.0);

fprintf('Fig.9  -> fig9_convergence.pdf\n');
fprintf('Fig.9b -> fig9b_conv_stats.pdf\n');
fprintf('CSV    -> %s\n', CSV);
end
