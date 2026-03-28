function run_fig7_robustness(fast, use_par)
%  use_par : logical (default false) -- enables parfor across MC trials
%RUN_FIG7_ROBUSTNESS  Fig. 7 -- Robustness: exact USW data, Fresnel estimators.
%  Results written to nf_simulation_results.csv after each SNR point.
if nargin < 1; fast = false; end
if nargin < 2; use_par = false; end

CSV = 'nf_simulation_results.csv';

P      = nf_params();
P.M    = 64;  P = nf_update_derived_pub(P);
P.N_RF = 8;   P.N = 64;   P.d = 3;
P.N_MC = 400; if fast; P.N_MC = 20; end

% Explicit 9-point sweep matching Fig.2
SNR_vec = [-15, -10:5:20, 25];
n_snr   = numel(SNR_vec);
methods = {'CL-KL','P-SOMP','DL-OMP','DFrFT-NOMP','BF-SOMP'};
n_meth  = numel(methods);
assert(n_meth==5,'expected 5 methods');
NMSE_db = nan(n_meth,n_snr);  NMSE_std = nan(n_meth,n_snr);

fprintf('=== Fig.7: Robustness (exact USW truth, Fresnel est.)  N_MC=%d ===\n', P.N_MC);

for si = 1:n_snr
    SNR = SNR_vec(si);
    fprintf('SNR = %+3d dB  ...', SNR);

    % Exact USW truth; all estimators use Fresnel (deliberate mismatch)
    [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR, true);

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
    fprintf(' done. CL-KL NMSE = %.1f dB\n', NMSE_db(1,si));

    save_results_csv(CSV,'Fig7','Robustness_USW_mismatch','mismatch_robustness', ...
        P.M,P.N_RF,P.N,P.d,SNR,'SNR_dB',SNR,'exact_USW', ...
        methods,nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,[],[],P.N_MC,P.r_RD,c_iters,c_conv,c_N0);
end


plot_fig(SNR_vec, NMSE_db, NMSE_std, methods, ...
    {'Fig.7: Robustness -- Exact USW Data, Fresnel Estimators'; ...
     '(deliberate truth-model mismatch)'}, ...
    'SNR (dB)', 'Fig7_Robustness', ...
    'NumColumns', 2, 'LegTitle', '');
xlim([min(SNR_vec)-1, max(SNR_vec)+1]);
% Single-column export (side-by-side with Fig.8 in LaTeX)
nf_export_fig(gcf, 'fig7_robustness', 'single', 'Height', 8.0);
fprintf('Fig.7 -> fig7_robustness.pdf  |  CSV -> %s\n', CSV);
end


% ====================================================================
function res = mc_trial(Xm, Hm, tht, rt, Wm, Ym, Rm, P)
%MC_TRIAL  Run CL-KL, P-SOMP, DL-OMP, DFrFT-NOMP (all Fresnel est.).
nm=nan(5,1); rth=nan(5,1); rr=nan(5,1); fl=nan(5,1);
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

try; [th,rh]=nf_dfrft_nomp(Xm,Wm,P);
    [nm(4),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(4)=a*180/pi; rr(4)=b; fl(4)=c*100;
catch; nm(4)=1; fl(4)=100; end


try; [th,rh]=nf_bfsomp(Xm,Wm,P);
    [nm(5),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
    rth(5)=a*180/pi; rr(5)=b; fl(5)=c*100;
catch; nm(5)=1; fl(5)=100; end
res = struct('nm',nm,'rth',rth,'rr',rr,'fl',fl, ...
    'ci_iter',ci_iter,'ci_conv',ci_conv,'ci_n0',ci_n0);
end
