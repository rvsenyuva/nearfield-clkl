function run_fig8_vard(fast, use_par)
%RUN_FIG8_VARD  Fig. 8 -- NMSE vs number of paths d.
%
%  New figure added in Step 5 of the TWC submission improvement plan.
%
%  PURPOSE:
%  Reviewers will ask "does this work for d != 3?" since every other figure
%  fixes d=3.  This figure sweeps d in {1,2,3,4,5} at fixed SNR=10dB and
%  shows all six estimators plus the CRB.
%
%  DESIGN NOTES:
%  -------------
%  1. Identifiability: compressed-domain methods (CL-KL, P-SOMP, BF-SOMP)
%     require N_RF >= 2d+1 for the N_RF x N_RF covariance to uniquely
%     identify d paths.  With N_RF=8: identifiable for d<=3, under-determined
%     for d>=4.  The figure exposes this degradation naturally -- no special
%     treatment needed.  A vertical dashed line at d=3 (N_RF/2 rounded down)
%     marks the identifiability boundary.
%
%  2. BF-SOMP adaptive: BF-SOMP estimates d from data (residual threshold).
%     We add a second BF-SOMP curve "BF-SOMP (adaptive)" where d_cap is set
%     to a loose upper bound (d_cap = 2*P.d) so BF-SOMP chooses its own
%     sparsity level.  This showcases its key advantage over the other methods
%     which all require d to be known exactly.
%
%  3. CRB: computed at each d value via nf_crb.m, averaged over MC trials.
%     May become ill-conditioned at d=4,5 with N_RF=8 -- regularisation in
%     nf_crb handles this gracefully.
%
%  OUTPUTS: fig8_vard.pdf, fig8b_vard_rmse.pdf, CSV entries tagged 'Fig8'.
%
%  P9.2 changes (19 April 2026):
%    - CRB loop now calls nf_crb(..., 'both') and extracts both
%      compressed and full-array CRB structs.
%    - save_results_csv now receives crb_theta_full_deg and crb_r_full_m
%      as new args 28-29 (30-column schema, P9.1).
%    - Fig.8b RMSE panels gain a second dash-dot CRB line for the
%      full-array bound; existing compressed CRB line style is unchanged.
%    - All method simulation code is bit-identical to Run 4.
if nargin < 1; fast = false; end
if nargin < 2; use_par = false; end
rng(42, 'twister');  % Fixed seed: ensures bit-exact reproducibility across runs

CSV    = 'nf_simulation_results.csv';
P_base = nf_params();
P_base.M    = 64;  P_base = nf_update_derived_pub(P_base);
P_base.N_RF = 8;   P_base.N = 64;
P_base.N_MC = 400; if fast; P_base.N_MC = 20; end
SNR_fix  = 10;

d_vec  = [1, 2, 3, 4, 5];
n_d    = numel(d_vec);
% d_cap for BF-SOMP adaptive: use 2x true d, minimum 3
d_cap_fac = 2;

methods = {'CL-KL','P-SOMP','DL-OMP','MUSIC+Tri','DFrFT-NOMP','BF-SOMP','BF-SOMP (adpt)'};
n_meth  = numel(methods);   % 7

NMSE_db   = nan(n_meth, n_d);
NMSE_std  = nan(n_meth, n_d);
RMSE_th   = nan(n_meth, n_d);
RMSE_r    = nan(n_meth, n_d);

% CRB arrays — compressed and full-array (median over MC trials)
CRB_theta      = nan(1, n_d);   % compressed CRB angle [deg]
CRB_r_m        = nan(1, n_d);   % compressed CRB range [m]
CRB_theta_full = nan(1, n_d);   % full-array CRB angle [deg]  (P9.2)
CRB_r_full     = nan(1, n_d);   % full-array CRB range [m]    (P9.2)

fprintf('=== Fig.8: NMSE vs d (N_MC=%d, SNR=%ddB) ===\n', P_base.N_MC, SNR_fix);
fprintf('    N_RF=%d: identifiable for d <= %d (N_RF >= 2d+1)\n', ...
    P_base.N_RF, floor((P_base.N_RF-1)/2));

for di = 1:n_d
    d_true = d_vec(di);
    P      = P_base;
    P.d    = d_true;
    P      = nf_update_derived_pub(P);
    d_cap  = max(d_true, min(d_cap_fac*d_true, P.N_RF-1));

    fprintf('d = %d  (d_cap_adaptive=%d)  ...', d_true, d_cap);

    [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR_fix, true);

    % -----------------------------------------------------------------
    % CRB: call nf_crb with mode='both' to obtain compressed AND
    % full-array bounds in one pass.  50 trials sufficient.
    % -----------------------------------------------------------------
    crb_th_mc      = nan(1, P.N_MC);
    crb_r_mc       = nan(1, P.N_MC);
    crb_th_full_mc = nan(1, P.N_MC);
    crb_r_full_mc  = nan(1, P.N_MC);

    for mc_crb = 1:min(P.N_MC, 50)   % 50 trials sufficient for CRB estimate
        [X_crb,~,~,th_crb,r_crb,p_crb,N0_crb] = nf_gen_channel(P, SNR_fix, true);
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
    CRB_theta(di)      = nanmedian(crb_th_mc);
    CRB_r_m(di)        = nanmedian(crb_r_mc);
    CRB_theta_full(di) = nanmedian(crb_th_full_mc);
    CRB_r_full(di)     = nanmedian(crb_r_full_mc);

    nmse_mc    = nan(n_meth, P.N_MC);
    rmse_th_mc = nan(n_meth, P.N_MC);
    rmse_r_mc  = nan(n_meth, P.N_MC);
    fail_mc    = nan(n_meth, P.N_MC);
    clkl_iters_mc = nan(1, P.N_MC);
    clkl_conv_mc  = nan(1, P.N_MC);
    clkl_N0_mc    = nan(1, P.N_MC);

    res_all = cell(P.N_MC, 1);
    if use_par
        parfor mc = 1:P.N_MC
            res_all{mc} = mc_trial(X_all{mc},H_all{mc},th_all{mc},r_all{mc}, ...
                W_all{mc},Y_all{mc},R_all{mc},P,d_cap);
        end
    else
        for mc = 1:P.N_MC
            res_all{mc} = mc_trial(X_all{mc},H_all{mc},th_all{mc},r_all{mc}, ...
                W_all{mc},Y_all{mc},R_all{mc},P,d_cap);
        end
    end
    for mc = 1:P.N_MC
        r = res_all{mc};
        nmse_mc(:,mc)    = r.nm;  rmse_th_mc(:,mc) = r.rth;
        rmse_r_mc(:,mc)  = r.rr;  fail_mc(:,mc)    = r.fl;
        clkl_iters_mc(mc)=r.ci_iter; clkl_conv_mc(mc)=r.ci_conv;
        clkl_N0_mc(mc)   =r.ci_n0;
    end

    [nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,c_iters,c_conv,c_N0] = ...
        aggregate(nmse_mc,rmse_th_mc,rmse_r_mc,fail_mc,n_meth, ...
                  clkl_iters_mc,clkl_conv_mc,clkl_N0_mc);
    NMSE_db(:,di)  = nmse_pt;
    NMSE_std(:,di) = nstd_pt;
    RMSE_th(:,di)  = rth_pt;
    RMSE_r(:,di)   = rr_pt;
    fprintf(' done. CL-KL NMSE=%.1fdB, CRB_th=%.3fdeg, CRB_r=%.2fm\n', ...
        nmse_pt(1), CRB_theta(di), CRB_r_m(di));

    % -----------------------------------------------------------------
    % save_results_csv: pass both compressed and full-array CRBs
    % (new args 28-29 per 30-column schema from P9.1)
    % -----------------------------------------------------------------
    save_results_csv(CSV,'Fig8','NMSE_vs_d','vard_sweep', ...
        P.M,P.N_RF,P.N,d_true,SNR_fix,'d',d_true,'exact_USW', ...
        methods,nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,[],[],P.N_MC,P.r_RD, ...
        c_iters,c_conv,c_N0, CRB_theta(di), CRB_r_m(di), ...
        CRB_theta_full(di), CRB_r_full(di));
end


% ---- Plot -----------------------------------------------------------
styles7 = {'-o','-s','--^','-.d','-v','-.h','--p'};
colors7  = {[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19], ...
            [0.49 0.18 0.56],[0.93 0.69 0.13],[0.30 0.57 0.43], ...
            [0.30 0.57 0.43]};

d_id = floor((P_base.N_RF-1)/2);

% --- Fig.8: NMSE vs d (single panel, legend southoutside) ------------
fig8 = figure('Name','Fig8_NMSE_vs_d','Position',[100 100 560 420]);
for mi = 1:n_meth-1
    errorbar(d_vec, NMSE_db(mi,:), NMSE_std(mi,:)/2, ...
        styles7{mi},'Color',colors7{mi},'LineWidth',1.5,'MarkerSize',7, ...
        'DisplayName',methods{mi});
    hold on;
end
% BF-SOMP adaptive: dashed, same colour
plot(d_vec, NMSE_db(n_meth,:), styles7{n_meth}, 'Color', colors7{n_meth}, ...
    'LineWidth',1.5,'MarkerSize',7,'LineStyle','--','DisplayName',methods{n_meth});
xline(d_id+0.5,'--k','LineWidth',1.2, ...
    'Label','N_{RF}/2 ident. limit','LabelVerticalAlignment','bottom', ...
    'HandleVisibility','off');
xlabel('Number of paths d','FontSize',12);
ylabel('NMSE (dB)','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on; set(gca,'XTick',d_vec,'FontSize',11);
% NumColumns=2 -> 4 rows; narrower box fits single-column width
leg8 = legend('Orientation','horizontal','NumColumns',2, ...
              'Location','southoutside','FontSize',8);
leg8.Title.String = '';
% Height=9.5 cm: sufficient for 4-row legend below axes
nf_export_fig(fig8, 'fig8_vard', 'single', 'Height', 9.5);

% --- Fig.8b: RMSE vs d with dual CRB lines (two subplots, unified legend) ---
%
%  P9.2: two CRB lines per panel
%    - compressed CRB : black dashed ('k--'), label 'CRB (compressed)'
%      (style unchanged from Run 4 — was labelled 'CRB')
%    - full-array CRB : grey dash-dot (colour [0.5 0.5 0.5], '-.',
%      LineWidth 1.5), label 'CRB (full)'
%
fig8b = figure('Name','Fig8b_RMSE_vs_d','Position',[120 120 900 380]);

subplot(1,2,1);
for mi = 1:n_meth
    plot(d_vec, RMSE_th(mi,:), styles7{mi},'Color',colors7{mi}, ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi}); hold on;
end
% Compressed CRB — black dashed (unchanged from Run 4, label updated)
plot(d_vec, CRB_theta,'k--','LineWidth',2,'DisplayName','CRB (compressed)');
% Full-array CRB — grey dash-dot (P9.2 addition)
plot(d_vec, CRB_theta_full,'-.','Color',[0.5 0.5 0.5],'LineWidth',1.5, ...
    'DisplayName','CRB (full)');
xline(d_id+0.5,'--','Color',[0.5 0.5 0.5],'LineWidth',1.0,'HandleVisibility','off');
xlabel('d','FontSize',12); ylabel('RMSE(\theta) [deg]','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on; set(gca,'XTick',d_vec);

subplot(1,2,2);
for mi = 1:n_meth
    semilogy(d_vec, RMSE_r(mi,:), styles7{mi},'Color',colors7{mi}, ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi}); hold on;
end
% Compressed CRB — black dashed (unchanged from Run 4, label updated)
semilogy(d_vec, CRB_r_m,'k--','LineWidth',2,'DisplayName','CRB (compressed)');
% Full-array CRB — grey dash-dot (P9.2 addition)
semilogy(d_vec, CRB_r_full,'-.','Color',[0.5 0.5 0.5],'LineWidth',1.5, ...
    'DisplayName','CRB (full)');
xline(d_id+0.5,'--','Color',[0.5 0.5 0.5],'LineWidth',1.0,'HandleVisibility','off');
xlabel('d','FontSize',12); ylabel('RMSE(r) [m]','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on; set(gca,'XTick',d_vec);

% sgtitle removed (P6.1): LaTeX caption is canonical.
% Unified legend: 7 methods + 2 CRB lines = 9 entries
methods8b = [methods, {'CRB (compressed)', 'CRB (full)'}];
styles8b  = [styles7, {'--', '-.'}];
colors8b  = [colors7, {[0 0 0], [0.5 0.5 0.5]}];
nf_add_legend(fig8b, methods8b, styles8b, colors8b, ...
    'FontSize', 8, 'filled_idx', [1 2], 'NumColumns', 4);
% Height=9.5 cm for 2-3 row legend (9 entries at NumColumns=4)
nf_export_fig(fig8b, 'fig8b_vard_rmse', 'double', 'Height', 9.5);

fprintf('Fig.8  -> fig8_vard.pdf\n');
fprintf('Fig.8b -> fig8b_vard_rmse.pdf\n');
fprintf('CSV    -> %s\n', CSV);
end

% ====================================================================
function res = mc_trial(Xm, Hm, tht, rt, Wm, Ym, Rm, P, d_cap)
%MC_TRIAL  Run all 7 estimators (6 standard + BF-SOMP adaptive).
%  d_cap: loose sparsity cap for the adaptive BF-SOMP variant.
nm=nan(7,1); rth=nan(7,1); rr=nan(7,1); fl=nan(7,1);
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

% BF-SOMP with EXACT d (same cap as all others)
try; [th,rh]=nf_bfsomp(Xm,Wm,P);
    [nm(6),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(6)=a*180/pi; rr(6)=b; fl(6)=c*100;
catch; nm(6)=1; fl(6)=100; end

% BF-SOMP ADAPTIVE: loose d_cap = 2*d_true (estimates d from residual)
P_adpt   = P;
P_adpt.d = d_cap;
try; [th,rh]=nf_bfsomp(Xm,Wm,P_adpt);
    % Evaluate with TRUE d for fair NMSE computation
    [nm(7),a,b,c]=nf_metrics(th(1:P.d),rh(1:P.d),Ym,Wm,Hm,tht,rt,P);
    rth(7)=a*180/pi; rr(7)=b; fl(7)=c*100;
catch; nm(7)=1; fl(7)=100; end

res = struct('nm',nm,'rth',rth,'rr',rr,'fl',fl, ...
    'ci_iter',ci_iter,'ci_conv',ci_conv,'ci_n0',ci_n0);
end
