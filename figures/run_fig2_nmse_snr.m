function run_fig2_nmse_snr(fast, use_par)
%  use_par : logical  (default false) -- pass true to use parfor across MC trials
%RUN_FIG2_NMSE_SNR  Fig. 2 -- NMSE vs SNR (primary figure).
%  Results written to nf_simulation_results.csv after each SNR point.
%
%  P9.2 changes (19 April 2026):
%    - CRB loop now calls nf_crb(..., 'both') and extracts both
%      compressed and full-array CRB structs.
%    - save_results_csv now receives crb_theta_full_deg and crb_r_full_m
%      as new args 28-29 (30-column schema, P9.1).
%    - Fig.2b RMSE panels gain a second dash-dot CRB line for the
%      full-array bound; existing compressed CRB line style is unchanged.
%    - All method simulation code is bit-identical to Run 4.
if nargin < 1; fast = false; end
if nargin < 2; use_par = false; end
rng(42, 'twister');  % Fixed seed: ensures bit-exact reproducibility across runs

CSV = 'nf_simulation_results.csv';

P      = nf_params();
P.M    = 64;   P = nf_update_derived_pub(P);
P.N_RF = 8;    P.N = 64;    P.d = 3;
P.N_MC = 400;  if fast; P.N_MC = 20; end

% Explicit 9-point sweep: [-15,-10,-5,0,+5,+10,+15,+20,+25] dB
% (override P.SNR_dB_vec to ensure -15 and +25 are included)
SNR_vec = [-15, -10:5:20, 25];
n_snr   = numel(SNR_vec);
methods = {'CL-KL','P-SOMP','DL-OMP','MUSIC+Tri','DFrFT-NOMP','BF-SOMP'};
n_meth  = numel(methods);
assert(n_meth==6,'expected 6 methods');
NMSE_db = nan(n_meth, n_snr);  NMSE_std = nan(n_meth, n_snr);

% CRB arrays — compressed and full-array (median over MC trials)
CRB_theta_deg      = nan(1, n_snr);   % compressed CRB angle [deg]
CRB_r_m            = nan(1, n_snr);   % compressed CRB range [m]
CRB_theta_full_deg = nan(1, n_snr);   % full-array CRB angle [deg]  (P9.2)
CRB_r_full_m       = nan(1, n_snr);   % full-array CRB range [m]    (P9.2)

RMSE_th_all   = nan(n_meth, n_snr);   % for Fig.2b
RMSE_r_all    = nan(n_meth, n_snr);   % for Fig.2b

fprintf('=== Fig.2: NMSE vs SNR  (N_MC=%d) ===\n', P.N_MC);

for si = 1:n_snr
    SNR = SNR_vec(si);
    fprintf('SNR = %+3d dB  ... ', SNR);

    [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR, true);

    % -----------------------------------------------------------------
    % CRB: call nf_crb with mode='both' to obtain compressed AND
    % full-array bounds in one pass.  Average over MC trials.
    % -----------------------------------------------------------------
    crb_th_mc      = nan(1, P.N_MC);
    crb_r_mc       = nan(1, P.N_MC);
    crb_th_full_mc = nan(1, P.N_MC);
    crb_r_full_mc  = nan(1, P.N_MC);

    for mc_crb = 1:P.N_MC
        [X_crb,~,~,th_crb,r_crb,p_crb,N0_crb] = ...
            nf_gen_channel(P, SNR, true);
        [W_crb,~,~] = nf_hybrid_combiner(X_crb, P);
        try
            % 'both' mode returns a struct with four fields
            crb_s = nf_crb(th_crb, r_crb, p_crb, N0_crb, W_crb, P, 'both');
            crb_th_mc(mc_crb)      = crb_s.theta_comp;
            crb_r_mc(mc_crb)       = crb_s.r_comp;
            crb_th_full_mc(mc_crb) = crb_s.theta_full;
            crb_r_full_mc(mc_crb)  = crb_s.r_full;
        catch; end
    end
    % Median: robust to ill-conditioned outliers (same aggregation as Run 4)
    CRB_theta_deg(si)      = nanmedian(crb_th_mc);
    CRB_r_m(si)            = nanmedian(crb_r_mc);
    CRB_theta_full_deg(si) = nanmedian(crb_th_full_mc);
    CRB_r_full_m(si)       = nanmedian(crb_r_full_mc);

    nmse_mc    = nan(n_meth,P.N_MC);
    rmse_th_mc = nan(n_meth,P.N_MC);
    rmse_r_mc  = nan(n_meth,P.N_MC);
    fail_mc    = nan(n_meth,P.N_MC);
    clkl_iters_mc = nan(1,P.N_MC);
    clkl_conv_mc  = nan(1,P.N_MC);
    clkl_N0_mc    = nan(1,P.N_MC);

    res_all = cell(P.N_MC, 1);
    if use_par
        parfor mc = 1:P.N_MC
            res_all{mc} = mc_trial(X_all{mc},H_all{mc},th_all{mc},r_all{mc},W_all{mc},Y_all{mc},R_all{mc},P);
        end
    else
        for mc = 1:P.N_MC
            res_all{mc} = mc_trial(X_all{mc},H_all{mc},th_all{mc},r_all{mc},W_all{mc},Y_all{mc},R_all{mc},P);
        end
    end
    for mc = 1:P.N_MC
        r = res_all{mc};
        nmse_mc(:,mc)=r.nm; rmse_th_mc(:,mc)=r.rth;
        rmse_r_mc(:,mc)=r.rr; fail_mc(:,mc)=r.fl;
        clkl_iters_mc(mc)=r.ci_iter; clkl_conv_mc(mc)=r.ci_conv; clkl_N0_mc(mc)=r.ci_n0;
    end

    [nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,c_iters,c_conv,c_N0] = aggregate(nmse_mc,rmse_th_mc,rmse_r_mc,fail_mc,n_meth,clkl_iters_mc,clkl_conv_mc,clkl_N0_mc);
    NMSE_db(:,si)=nmse_pt;  NMSE_std(:,si)=nstd_pt;
    RMSE_th_all(:,si)=rth_pt;  RMSE_r_all(:,si)=rr_pt;
    fprintf('done. CL-KL NMSE = %.1f dB\n', NMSE_db(1,si));

    % -----------------------------------------------------------------
    % save_results_csv: pass both compressed and full-array CRBs
    % (new args 28-29 per 30-column schema from P9.1)
    % -----------------------------------------------------------------
    save_results_csv(CSV,'Fig2','NMSE_vs_SNR','primary', ...
        P.M,P.N_RF,P.N,P.d,SNR,'SNR_dB',SNR,'exact_USW', ...
        methods,nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,[],[],P.N_MC,P.r_RD, ...
        c_iters,c_conv,c_N0, CRB_theta_deg(si), CRB_r_m(si), ...
        CRB_theta_full_deg(si), CRB_r_full_m(si));
end


plot_fig(SNR_vec, NMSE_db, NMSE_std, methods, ...
    'Fig.2: Channel NMSE vs SNR', 'SNR (dB)', 'Fig2_NMSE_vs_SNR');

% --- R1: Annotate CL-KL +25 dB OMP warm-start limitation ---
% Find the +25 dB point index
idx25 = find(SNR_vec == 25, 1);
if ~isempty(idx25)
    ax = gca;
    hold(ax, 'on');
    % Plot hollow circle marker on the CL-KL curve at +25 dB
    plot(ax, SNR_vec(idx25), NMSE_db(1,idx25), 'ko', ...
        'MarkerSize', 10, 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'none', 'HandleVisibility', 'off');
    % Text annotation: use latex interpreter so \dagger renders correctly
    text(ax, SNR_vec(idx25) - 1.2, NMSE_db(1,idx25) + 2.0, ...
        '$\dagger\dagger$ OMP limitation', ...
        'FontSize', 7, 'HorizontalAlignment', 'right', ...
        'Interpreter', 'latex');
end

% Height 9.5 cm: +2 cm vs original 7.5 to accommodate below-axes legend
nf_export_fig(gcf, 'fig2_nmse_snr', 'double', 'Height', 9.5);

% ---- Fig.2b: RMSE_theta and RMSE_r with dual CRB lines and unified legend ---
%
%  P9.2: two CRB lines per panel
%    - compressed CRB : black dashed ('k--'), label 'CRB (compressed)'
%      (style unchanged from Run 4 — was labelled 'CRB')
%    - full-array CRB : grey dash-dot (colour [0.5 0.5 0.5], '--.',
%      LineWidth 1.5), label 'CRB (full)'
%
styles6={'-o','-s','--^','-.d','-v','-.h'};
colors6={[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19], ...
         [0.49 0.18 0.56],[0.93 0.69 0.13],[0.30 0.57 0.43]};
fig2b = figure('Name','Fig2b_RMSE_CRB','Position',[120 120 900 380]);
RMSE_th_db = RMSE_th_all;  RMSE_r_db = RMSE_r_all;

subplot(1,2,1);
for mi=1:n_meth
    plot(SNR_vec,RMSE_th_db(mi,:),styles6{mi},'Color',colors6{mi}, ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi}); hold on;
end
% Compressed CRB — black dashed (unchanged from Run 4, label updated)
plot(SNR_vec,CRB_theta_deg,'k--','LineWidth',2,'DisplayName','CRB (compressed)');
% Full-array CRB — grey dash-dot (P9.2 addition)
plot(SNR_vec,CRB_theta_full_deg,'--','Color',[0.5 0.5 0.5],'LineWidth',1.5, ...
    'LineStyle','-.','DisplayName','CRB (full)');
xlabel('SNR (dB)','FontSize',12); ylabel('RMSE(\theta) [deg]','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on;

subplot(1,2,2);
for mi=1:n_meth
    semilogy(SNR_vec,RMSE_r_db(mi,:),styles6{mi},'Color',colors6{mi}, ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi}); hold on;
end
% Compressed CRB — black dashed (unchanged from Run 4, label updated)
semilogy(SNR_vec,CRB_r_m,'k--','LineWidth',2,'DisplayName','CRB (compressed)');
% Full-array CRB — grey dash-dot (P9.2 addition)
semilogy(SNR_vec,CRB_r_full_m,'-.','Color',[0.5 0.5 0.5],'LineWidth',1.5, ...
    'DisplayName','CRB (full)');
xlabel('SNR (dB)','FontSize',12); ylabel('RMSE(r) [m]','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on;

% sgtitle removed (P6.1): LaTeX caption is canonical.
% Unified legend: 8 entries (6 methods + 2 CRB lines)
methods_crb = [methods, {'CRB (compressed)', 'CRB (full)'}];
styles_crb  = [styles6, {'--', '-.'}];
colors_crb  = [colors6, {[0 0 0], [0.5 0.5 0.5]}];
nf_add_legend(fig2b, methods_crb, styles_crb, colors_crb, ...
    'FontSize', 8, 'filled_idx', [1 2], 'NumColumns', 4);
% Height 8.5 cm: +1.5 cm vs original 7.0 to accommodate unified legend
nf_export_fig(fig2b, 'fig2b_rmse_crb', 'double', 'Height', 8.5);
fprintf('Fig.2b -> fig2b_rmse_crb.pdf\n');
fprintf('Fig.2 -> fig2_nmse_snr.pdf  |  CSV -> %s\n', CSV);
end

% ====================================================================
function res = mc_trial(Xm, Hm, tht, rt, Wm, Ym, Rm, P)
%MC_TRIAL  Run all five estimators on one pre-generated trial.
nm=nan(6,1); rth=nan(6,1); rr=nan(6,1); fl=nan(6,1);
ci_iter=nan; ci_conv=nan; ci_n0=nan;

try; [th,rh,~,~,ci]=nf_clkl(Rm,Wm,P);
    ci_iter=ci.n_iter; ci_conv=ci.converged; ci_n0=ci.N0_hat;
    [nm(1),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(1)=a*180/pi; rr(1)=b; fl(1)=c*100;
catch; nm(1)=1; fl(1)=100; end

try; [th,rh]=nf_psomp(Rm,Wm,P);
    [nm(2),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(2)=a*180/pi; rr(2)=b; fl(2)=c*100;
catch; nm(2)=1; fl(2)=100; end

try; [th,rh]=nf_zhang(Xm,Wm,P);
    [nm(3),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(3)=a*180/pi; rr(3)=b; fl(3)=c*100;
catch; nm(3)=1; fl(3)=100; end

try; [th,rh]=nf_music_tri(Xm,P);
    [nm(4),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(4)=a*180/pi; rr(4)=b; fl(4)=c*100;
catch; nm(4)=1; fl(4)=100; end

try; [th,rh]=nf_dfrft_nomp(Xm,Wm,P);
    [nm(5),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(5)=a*180/pi; rr(5)=b; fl(5)=c*100;
catch; nm(5)=1; fl(5)=100; end


try; [th,rh]=nf_bfsomp(Xm,Wm,P);
    [nm(6),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(6)=a*180/pi; rr(6)=b; fl(6)=c*100;
catch; nm(6)=1; fl(6)=100; end
res = struct('nm',nm,'rth',rth,'rr',rr,'fl',fl, ...
    'ci_iter',ci_iter,'ci_conv',ci_conv,'ci_n0',ci_n0);
end
