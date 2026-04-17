function run_fig6_runtime(fast, use_par)  %#ok<INUSD>
%  use_par is accepted but ignored -- Fig.6 always runs serially.
%  Timing measurements must reflect single-core cost; parfor would
%  distort the runtime vs M comparison that is the purpose of this figure.
%RUN_FIG6_RUNTIME  Fig. 6 -- Runtime vs M (computational scaling).
%  Results (runtime_s AND NMSE) written to nf_simulation_results.csv.
%  Serial MC intentional -- measures true single-core per-method cost.
if nargin < 1; fast = false; end
% use_par intentionally ignored -- serial timing is required for Fig.6
rng(42, 'twister');  % Fixed seed: ensures bit-exact reproducibility across runs

CSV = 'nf_simulation_results.csv';

P      = nf_params();
P.N_RF = 8;   P.N = 64;   P.d = 3;
P.N_MC = 50;  if fast; P.N_MC = 10; end
SNR_fix = 10;

M_vec   = [32, 64, 128, 256];   % M=256 added to demonstrate M-invariant CL-KL cost
n_M     = numel(M_vec);
methods = {'CL-KL','P-SOMP','DL-OMP','DFrFT-NOMP','BF-SOMP'};
n_meth  = numel(methods);
assert(n_meth==5,'expected 5 methods');
RT_mean = nan(n_meth,n_M);
NMSE_db = nan(n_meth,n_M);

fprintf('=== Fig.6: Runtime vs M  (N_MC=%d per point) ===\n', P.N_MC);

for mi_idx = 1:n_M
    P.M = M_vec(mi_idx);
    P   = nf_update_derived_pub(P);
    fprintf('M = %3d  ...', P.M);

    [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR_fix, true);

    rt_mc      = nan(n_meth,P.N_MC);
    nmse_mc    = nan(n_meth,P.N_MC);
    rmse_th_mc = nan(n_meth,P.N_MC);
    rmse_r_mc  = nan(n_meth,P.N_MC);
    fail_mc    = nan(n_meth,P.N_MC);
    % CL-KL convergence diagnostics (method 1 only)
    clkl_iters_mc = nan(1,P.N_MC);
    clkl_conv_mc  = nan(1,P.N_MC);
    clkl_N0_mc    = nan(1,P.N_MC);

    % Serial -- fair single-core timing comparison
    for mc = 1:P.N_MC
        Xm=X_all{mc}; Hm=H_all{mc}; tht=th_all{mc}; rt=r_all{mc};
        Wm=W_all{mc}; Ym=Y_all{mc}; Rm=R_all{mc};
        for mi=1:n_meth
            t0=tic;
            try
                switch mi
                    case 1; [th,rh,~,~,ci]=nf_clkl(Rm,Wm,P); clkl_iters_mc(mc)=ci.n_iter; clkl_conv_mc(mc)=ci.converged; clkl_N0_mc(mc)=ci.N0_hat;
                    case 2; [th,rh]=nf_psomp(Rm,Wm,P);
                    case 3; [th,rh]=nf_zhang(Xm,Wm,P);
                    case 4; [th,rh]=nf_dfrft_nomp(Xm,Wm,P);
                    case 5; [th,rh]=nf_bfsomp(Xm,Wm,P);  % BF-SOMP full-array
                end
                rt_mc(mi,mc)=toc(t0);
                [nmse_mc(mi,mc),a,b,c]=nf_metrics(th,rh,Ym,Wm,Hm,tht,rt,P);
                rmse_th_mc(mi,mc)=a*180/pi; rmse_r_mc(mi,mc)=b; fail_mc(mi,mc)=c*100;
            catch
                rt_mc(mi,mc)=toc(t0); nmse_mc(mi,mc)=1; fail_mc(mi,mc)=100;
            end
        end
    end

    [nmse_pt,~,rth_pt,rr_pt,fail_pt,c_iters,c_conv,c_N0] = aggregate(nmse_mc,rmse_th_mc,rmse_r_mc,fail_mc,n_meth,clkl_iters_mc,clkl_conv_mc,clkl_N0_mc);
    RT_mean(:,mi_idx)  = mean(rt_mc,2,'omitnan');
    NMSE_db(:,mi_idx)  = nmse_pt;
    fprintf(' done. CL-KL t=%.2fs, P-SOMP t=%.2fs, BF-SOMP t=%.2fs\n',RT_mean(1,mi_idx),RT_mean(2,mi_idx),RT_mean(5,mi_idx));

    % Write runtime in the runtime_s column; fixed_SNR_dB = SNR_fix
    save_results_csv(CSV,'Fig6','Runtime_vs_M','scaling', ...
        M_vec(mi_idx),P.N_RF,P.N,P.d,SNR_fix,'M',M_vec(mi_idx),'exact_USW', ...
        methods,nmse_pt,nan(n_meth,1),rth_pt,rr_pt,fail_pt,RT_mean(:,mi_idx),[],P.N_MC,P.r_RD,c_iters,c_conv,c_N0);
end


fig6 = figure('Name','Fig6_Runtime_vs_M','Position',[100 100 800 380]);
styles={'-o','-s','--^',':d','-.h'};
colors={[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19],[0.49 0.18 0.56],[0.30 0.57 0.43]};

subplot(1,2,1);
for mi=1:n_meth
    semilogy(M_vec,RT_mean(mi,:),styles{mi},'Color',colors{mi},...
        'LineWidth',1.5,'MarkerSize',8,'DisplayName',methods{mi}); hold on;
end
xlabel('Array size M','FontSize',12); ylabel('Runtime per trial [s]','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on;
xticks(M_vec);

subplot(1,2,2);
for mi=1:n_meth
    plot(M_vec,NMSE_db(mi,:),styles{mi},'Color',colors{mi},...
        'LineWidth',1.5,'MarkerSize',8,'DisplayName',methods{mi}); hold on;
end
xlabel('Array size M','FontSize',12); ylabel('NMSE (dB)','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on; xticks(M_vec);

% sgtitle removed (P6.1): LaTeX caption is canonical.
% Unified single legend: filled_idx=[1 2] for CL-KL, P-SOMP
nf_add_legend(fig6, methods, styles, colors, 'FontSize',8, 'filled_idx',[1 2]);
% Height 8.5 cm (+1.5 cm vs original 7.0)
nf_export_fig(gcf, 'fig6_runtime', 'double', 'Height', 8.5);
fprintf('Fig.6 -> fig6_runtime.pdf  |  CSV -> %s\n', CSV);
end
