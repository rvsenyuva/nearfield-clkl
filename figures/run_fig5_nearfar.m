function run_fig5_nearfar(fast, use_par)
%  use_par : logical (default false) -- enables parfor across MC trials
%RUN_FIG5_NEARFAR  Fig. 5 -- Near->Far transition sweep (Step 1 updated).
%
%  CHANGES in this version (Step 1 -- Physical Regime Audit):
%  -----------------------------------------------------------
%  1. rmax_fac_vec extended from 8 to 10 points, starting at 0.1*r_RD.
%     The new points [0.1, 0.2] cover the strong near-field regime inside
%     the EBRD (r_EBRD ~ 0.058*r_RD angle-averaged), where polar-grid
%     coherence is most severe and CL-KL curvature learning is most needed.
%
%  2. r_lo_fac held FIXED at P.r_lo_fac (0.05) across all sweep points.
%     Previously it was set proportionally to rmax_fac, which moved the
%     entire window and prevented testing the strong NF at large r_max.
%     Now r_lo_fac=0.05 always; only r_hi_fac=rmax_fac varies.
%     Sources are always drawn from [0.05, rmax_fac]*r_RD.
%
%  3. Vertical lines added to all three subplots:
%     - r_EBRD (angle-averaged) = 0.058*r_RD (Hussain et al. TWC 2025 Thm.2)
%       marks the onset of strong polar-domain sparsity.
%     - r_RD = 1.0*r_RD marks the classical Rayleigh distance.
%
%  Results written to nf_simulation_results.csv after each sweep point.
if nargin < 1; fast = false; end
if nargin < 2; use_par = false; end
rng(42, 'twister');  % Fixed seed: ensures bit-exact reproducibility across runs

CSV = 'nf_simulation_results.csv';

P      = nf_params();
P.M    = 64;  P = nf_update_derived_pub(P);
P.N_RF = 8;   P.N = 64;   P.d = 3;
P.N_MC = 400; if fast; P.N_MC = 20; end
SNR_fix = 10;

% Step 1: extended range sweep (10 points, starts inside EBRD at 0.1*r_RD)
rmax_fac_vec = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0];
r_lo_fac_fix = P.r_lo_fac;   % 0.05 -- fixed across all sweep points
n_r    = numel(rmax_fac_vec);
methods = {'CL-KL','P-SOMP','DL-OMP','DFrFT-NOMP','BF-SOMP'};
n_meth  = numel(methods);
assert(n_meth==5,'expected 5 methods');

RMSE_r  = nan(n_meth,n_r);
RMSE_th = nan(n_meth,n_r);
NMSE_db = nan(n_meth,n_r);
r_RD_base = P.r_RD;

fprintf('=== Fig.5: Near->Far transition  (N_MC=%d) ===\n', P.N_MC);
fprintf('    r_lo fixed = %.2f*r_RD = %.2f m\n', r_lo_fac_fix, r_lo_fac_fix*r_RD_base);
fprintf('    EBRD (angle-avg) ~ %.3f*r_RD = %.2f m\n', P.r_EBRD_fac, P.r_EBRD_fac*r_RD_base);

for ri = 1:n_r
    rmax_fac   = rmax_fac_vec(ri);
    % Step 1: r_lo_fac is FIXED; only r_hi_fac varies
    P.r_hi_fac = rmax_fac;
    P.r_lo_fac = r_lo_fac_fix;
    fprintf('r_max/r_RD = %.2f  [%.2fm, %.2fm]  ...', ...
        rmax_fac, P.r_lo_fac*r_RD_base, P.r_hi_fac*r_RD_base);

    [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR_fix, true);

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
            res_all{mc} = mc_trial(X_all{mc},H_all{mc},th_all{mc},r_all{mc},W_all{mc},Y_all{mc},R_all{mc},P,r_RD_base);
        end
    else
        for mc = 1:P.N_MC
            res_all{mc} = mc_trial(X_all{mc},H_all{mc},th_all{mc},r_all{mc},W_all{mc},Y_all{mc},R_all{mc},P,r_RD_base);
        end
    end
    for mc = 1:P.N_MC
        r = res_all{mc};
        nmse_mc(:,mc)=r.nm; rmse_th_mc(:,mc)=r.rth;
        rmse_r_mc(:,mc)=r.rr; fail_mc(:,mc)=r.fl;
        clkl_iters_mc(mc)=r.ci_iter; clkl_conv_mc(mc)=r.ci_conv; clkl_N0_mc(mc)=r.ci_n0;
    end

    [nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,c_iters,c_conv,c_N0] = aggregate(nmse_mc,rmse_th_mc,rmse_r_mc,fail_mc,n_meth,clkl_iters_mc,clkl_conv_mc,clkl_N0_mc);
    NMSE_db(:,ri)=nmse_pt;
    RMSE_th(:,ri)=rth_pt;
    RMSE_r(:,ri) =rr_pt;
    fprintf(' done. CL-KL RMSE_r=%.1fm, RMSE_th=%.2f deg\n', rr_pt(1), rth_pt(1));

    save_results_csv(CSV,'Fig5','NearFar_transition','robustness_sweep', ...
        P.M,P.N_RF,P.N,P.d,SNR_fix,'r_max_over_rRD',rmax_fac,'exact_USW', ...
        methods,nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,[],[],P.N_MC,P.r_RD,c_iters,c_conv,c_N0);
end


% ---- Plot -- unified legend below, no per-subplot legends ---------------
fig5 = figure('Name','Fig5_NearFar','Position',[100 100 900 380]);
styles={'-o','-s','--^',':d','-.h'};
colors={[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19],[0.49 0.18 0.56],[0.30 0.57 0.43]};
x_ax = rmax_fac_vec;
ebrd_fac = P.r_EBRD_fac;
rrd_fac  = 1.0;

subplot(1,3,1);
for mi=1:n_meth
    semilogy(x_ax,RMSE_r(mi,:),styles{mi},'Color',colors{mi},...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi}); hold on;
end
xlabel('r_{max}/r_{RD}','FontSize',11); ylabel('RMSE(r) [m]','FontSize',11);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on;
xline(ebrd_fac,'--','Color',[0.6 0.2 0.8],'LineWidth',1.2,...
    'Label','r_{EBRD}','LabelVerticalAlignment','bottom',...
    'LabelHorizontalAlignment','right','HandleVisibility','off');
xline(rrd_fac,'--k','LineWidth',1.2,'Label','r_{RD}',...
    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left',...
    'HandleVisibility','off');

subplot(1,3,2);
for mi=1:n_meth
    plot(x_ax,RMSE_th(mi,:),styles{mi},'Color',colors{mi},...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi}); hold on;
end
xlabel('r_{max}/r_{RD}','FontSize',11); ylabel('RMSE(\theta) [deg]','FontSize',11);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on;
xline(ebrd_fac,'--','Color',[0.6 0.2 0.8],'LineWidth',1.2,'HandleVisibility','off');
xline(rrd_fac,'--k','LineWidth',1.2,'HandleVisibility','off');

subplot(1,3,3);
for mi=1:n_meth
    plot(x_ax,NMSE_db(mi,:),styles{mi},'Color',colors{mi},...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi}); hold on;
end
xlabel('r_{max}/r_{RD}','FontSize',11); ylabel('NMSE (dB)','FontSize',11);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on;
xline(ebrd_fac,'--','Color',[0.6 0.2 0.8],'LineWidth',1.2,'Label','EBRD',...
    'LabelVerticalAlignment','top','LabelHorizontalAlignment','right',...
    'HandleVisibility','off');
xline(rrd_fac,'--k','LineWidth',1.2,'Label','r_{RD}',...
    'LabelVerticalAlignment','top','LabelHorizontalAlignment','left',...
    'HandleVisibility','off');

% sgtitle removed (P6.1): LaTeX caption is canonical.

% Unified single legend centred below all three panels
% filled_idx = [1 2]: CL-KL and P-SOMP are compressed-domain (filled markers)
nf_add_legend(fig5, methods, styles, colors, 'FontSize',8, 'filled_idx',[1 2]);
% Height 9.5 cm (+1.5 cm vs original 8.0) to accommodate unified legend
nf_export_fig(gcf, 'fig5_nearfar', 'double', 'Height', 9.0);
fprintf('Fig.5 -> fig5_nearfar.pdf  |  CSV -> %s\n', CSV);
end


% ====================================================================
function res = mc_trial(Xm, Hm, tht, rt, Wm, Ym, Rm, P, r_RD_base)
%MC_TRIAL  Run CL-KL, P-SOMP, DL-OMP, DFrFT-NOMP on one trial.
nm=nan(5,1); rth=nan(5,1); rr=nan(5,1); fl=nan(5,1);
ci_iter=nan; ci_conv=nan; ci_n0=nan;

try; [th,rh,~,~,ci]=nf_clkl(Rm,Wm,P);
    ci_iter=ci.n_iter; ci_conv=ci.converged; ci_n0=ci.N0_hat;
    [nm(1),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(1)=a*180/pi; rr(1)=b; fl(1)=c*100;
catch; nm(1)=1; rr(1)=r_RD_base; fl(1)=100; end

try; [th,rh]=nf_psomp(Rm,Wm,P);
    [nm(2),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(2)=a*180/pi; rr(2)=b; fl(2)=c*100;
catch; nm(2)=1; rr(2)=r_RD_base; fl(2)=100; end

try; [th,rh]=nf_zhang(Xm,Wm,P);
    [nm(3),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(3)=a*180/pi; rr(3)=b; fl(3)=c*100;
catch; nm(3)=1; rr(3)=r_RD_base; fl(3)=100; end

try; [th,rh]=nf_dfrft_nomp(Xm,Wm,P);
    [nm(4),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(4)=a*180/pi; rr(4)=b; fl(4)=c*100;
catch; nm(4)=1; rr(4)=r_RD_base; fl(4)=100; end


try; [th,rh]=nf_bfsomp(Xm,Wm,P);
    [nm(5),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(5)=a*180/pi; rr(5)=b; fl(5)=c*100;
catch; nm(5)=1; fl(5)=100; end
res = struct('nm',nm,'rth',rth,'rr',rr,'fl',fl, ...
    'ci_iter',ci_iter,'ci_conv',ci_conv,'ci_n0',ci_n0);
end
