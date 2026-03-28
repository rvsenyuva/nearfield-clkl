function run_fig4_nmse_n(fast, use_par)
%  use_par : logical  (default false) — pass true to use parfor across MC trials
%RUN_FIG4_NMSE_N  Fig. 4 — NMSE vs N (snapshot sweep).
%  Results written to nf_simulation_results.csv after each N point.
if nargin < 1; fast = false; end
if nargin < 2; use_par = false; end

CSV = 'nf_simulation_results.csv';

P      = nf_params();
P.M    = 64;  P = nf_update_derived_pub(P);
P.N_RF = 8;   P.d = 3;
P.N_MC = 200; if fast; P.N_MC = 20; end
SNR_fix = 10;

N_vec  = [16, 32, 64, 128];
n_n    = numel(N_vec);
methods = {'CL-KL','P-SOMP','DL-OMP','MUSIC+Tri','DFrFT-NOMP','BF-SOMP'};
n_meth  = numel(methods);
assert(n_meth==6,'expected 6 methods');
NMSE_db = nan(n_meth,n_n);  NMSE_std = nan(n_meth,n_n);

fprintf('=== Fig.4: NMSE vs N  (SNR=%ddB, N_MC=%d) ===\n', SNR_fix, P.N_MC);

for ni = 1:n_n
    P.N = N_vec(ni);
    fprintf('N = %3d  ... ', P.N);

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
    NMSE_db(:,ni)=nmse_pt;  NMSE_std(:,ni)=nstd_pt;
    fprintf('done. CL-KL NMSE = %.1f dB\n', NMSE_db(1,ni));

    save_results_csv(CSV,'Fig4','NMSE_vs_N','snapshot_sweep', ...
        P.M,P.N_RF,N_vec(ni),P.d,SNR_fix,'N',N_vec(ni),'exact_USW', ...
        methods,nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,[],[],P.N_MC,P.r_RD,c_iters,c_conv,c_N0);
end


plot_fig(N_vec, NMSE_db, NMSE_std, methods, ...
    'Fig.4: NMSE vs Number of Snapshots', 'Snapshots N', 'Fig4_NMSE_vs_N', ...
    'NumColumns', 2, 'LegTitle', '');
set(gca,'XScale','log','XTick',N_vec);
% Single-column export (placed side-by-side with Fig.3 in LaTeX)
nf_export_fig(gcf, 'fig4_nmse_n', 'single', 'Height', 8.5);
fprintf('Fig.4 -> fig4_nmse_n.pdf  |  CSV -> %s\n', CSV);
end


% ====================================================================
function res = mc_trial(Xm, Hm, tht, rt, Wm, Ym, Rm, P)
%MC_TRIAL  Run all six estimators on one pre-generated trial.
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

try; [th,rh]=nf_bfsomp(Xm,Wm,P);   % full-array input
    [nm(6),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(6)=a*180/pi; rr(6)=b; fl(6)=c*100;
catch; nm(6)=1; fl(6)=100; end

res = struct('nm',nm,'rth',rth,'rr',rr,'fl',fl, ...
    'ci_iter',ci_iter,'ci_conv',ci_conv,'ci_n0',ci_n0);
end
