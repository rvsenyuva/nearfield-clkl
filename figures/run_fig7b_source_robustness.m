function run_fig7b_source_robustness(fast, use_par)
%RUN_FIG7B_SOURCE_ROBUSTNESS  Fig. 7b -- Source model robustness.
%
%  New figure added in Step 7 of the TWC submission improvement plan.
%
%  PURPOSE:
%  Reviewers familiar with Wu et al. (GOES, DSP 2020) and Zuo et al.
%  (Gridless SDP, ICASSP) will note that several near-field methods
%  require specific source distributions:
%    - GOES/W-GOES:  requires sub-Gaussian sources (4th-order cumulants)
%    - Zuo et al.:   uses second-order statistics only (any distribution)
%    - CL-KL:        KL objective between sample and model covariances;
%                    technically assumes Gaussian but since it fits only
%                    second-order statistics, it is robust to non-Gaussian.
%
%  This figure demonstrates CL-KL robustness by comparing NMSE under:
%    (a) Gaussian pilots:  s_l(t) ~ CN(0, p_l)
%    (b) QPSK pilots:      s_l(t) = sqrt(p_l) * exp(j*pi/4*(2k-1))
%
%  KEY PHYSICAL INSIGHT:
%  The signal covariance R_h = A*diag(p)*A^H depends only on p_l = E[|s_l|^2],
%  not on the distribution of s_l.  Covariance-fitting methods (CL-KL,
%  P-SOMP, BF-SOMP) therefore see identical second-order statistics for
%  both pilot types and should produce identical NMSE.  Methods that
%  exploit individual snapshots (DL-OMP, MUSIC+Tri, DFrFT-NOMP) receive
%  a different snapshot matrix X and may exhibit small differences.
%
%  EXPECTED RESULT:
%  CL-KL, P-SOMP, BF-SOMP: zero gap between Gaussian and QPSK curves.
%  DL-OMP, DFrFT-NOMP: small gap (< 1 dB) since their snapshot processing
%  implicitly exploits Gaussianity for optimal thresholding.
%  MUSIC+Tri: similar -- eigenvector decomposition is not distribution-free.
%
if nargin < 1; fast = false; end
if nargin < 2; use_par = false; end
rng(42, 'twister');  % Fixed seed: ensures bit-exact reproducibility across runs

CSV   = 'nf_simulation_results.csv';
P     = nf_params();
P.M   = 64;  P = nf_update_derived_pub(P);
P.N_RF = 8;  P.N = 64;  P.d = 3;
P.N_MC = 400; if fast; P.N_MC = 20; end

% Explicit 9-point sweep matching Fig.2
SNR_vec = [-15, -10:5:20, 25];
n_snr   = numel(SNR_vec);

methods = {'CL-KL','P-SOMP','DL-OMP','MUSIC+Tri','DFrFT-NOMP','BF-SOMP'};
n_meth  = numel(methods);

% Run both pilot types
pilot_types = {'gaussian', 'qpsk'};
n_pt = numel(pilot_types);
NMSE_all  = nan(n_meth, n_snr, n_pt);
NMSE_std_all = nan(n_meth, n_snr, n_pt);

for pti = 1:n_pt
    P.pilot_type = pilot_types{pti};
    fprintf('=== Fig.7b: Source robustness -- pilot=%s (N_MC=%d) ===\n', ...
        P.pilot_type, P.N_MC);

    for si = 1:n_snr
        SNR = SNR_vec(si);
        fprintf('  SNR = %+3d dB  ...', SNR);

        [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR, true);

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
                    W_all{mc},Y_all{mc},R_all{mc},P);
            end
        else
            for mc = 1:P.N_MC
                res_all{mc} = mc_trial(X_all{mc},H_all{mc},th_all{mc},r_all{mc}, ...
                    W_all{mc},Y_all{mc},R_all{mc},P);
            end
        end
        for mc = 1:P.N_MC
            r = res_all{mc};
            nmse_mc(:,mc)    = r.nm;  rmse_th_mc(:,mc) = r.rth;
            rmse_r_mc(:,mc)  = r.rr;  fail_mc(:,mc)    = r.fl;
            clkl_iters_mc(mc)= r.ci_iter; clkl_conv_mc(mc)= r.ci_conv;
            clkl_N0_mc(mc)   = r.ci_n0;
        end

        [nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,c_iters,c_conv,c_N0] = ...
            aggregate(nmse_mc,rmse_th_mc,rmse_r_mc,fail_mc,n_meth, ...
                      clkl_iters_mc,clkl_conv_mc,clkl_N0_mc);
        NMSE_all(:,si,pti)     = nmse_pt;
        NMSE_std_all(:,si,pti) = nstd_pt;
        fprintf(' done. CL-KL NMSE = %.1f dB\n', nmse_pt(1));

        % Tag experiment type in CSV
        exp_type = sprintf('source_robustness_%s', P.pilot_type);
        save_results_csv(CSV,'Fig7b','Source_model_robustness',exp_type, ...
            P.M,P.N_RF,P.N,P.d,SNR,'SNR_dB',SNR,P.pilot_type, ...
            methods,nmse_pt,nstd_pt,rth_pt,rr_pt,fail_pt,[],[],P.N_MC,P.r_RD, ...
            c_iters,c_conv,c_N0);
    end
end


% ---- Plot -- unified legend below, no per-subplot legends ---------------
styles6 = {'-o','-s','--^','-.d','-v','-.h'};
colors6  = {[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19], ...
            [0.49 0.18 0.56],[0.93 0.69 0.13],[0.30 0.57 0.43]};
compressed_idx = [1, 2];

fig7b = figure('Name','Fig7b_SourceRobust','Position',[100 100 1050 520]);

% --- Left panel: absolute NMSE for both pilot types ----------------------
subplot(1,2,1);
for mi = 1:n_meth
    mfc7 = colors6{mi};
    if ~ismember(mi, compressed_idx); mfc7 = 'none'; end
    errorbar(SNR_vec, NMSE_all(mi,:,1), NMSE_std_all(mi,:,1)/2, ...
        styles6{mi},'Color',colors6{mi},'MarkerFaceColor',mfc7, ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi});
    hold on;
    % QPSK: same colour, dashed, no errorbars, no legend entry
    plot(SNR_vec, NMSE_all(mi,:,2), '--','Color',colors6{mi}, ...
        'LineWidth',1.0,'HandleVisibility','off');
end
xlabel('SNR (dB)','FontSize',12); ylabel('NMSE (dB)','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on;

% --- Right panel: NMSE gap per method ------------------------------------
subplot(1,2,2);
gap_mat = NMSE_all(:,:,2) - NMSE_all(:,:,1);
for mi = 1:n_meth
    plot(SNR_vec, gap_mat(mi,:), styles6{mi},'Color',colors6{mi}, ...
        'LineWidth',1.5,'MarkerSize',7,'DisplayName',methods{mi});
    hold on;
end
yline(0,'--k','LineWidth',1.2,'HandleVisibility','off');
xlabel('SNR (dB)','FontSize',12);
ylabel('NMSE gap [dB] (QPSK - Gauss)','FontSize',12);
% title() call removed (P6.1): LaTeX caption is canonical.
grid on; ylim([-3 3]);

% sgtitle removed (P6.1): LaTeX caption is canonical.

% Unified single legend centred below both panels (methods only; QPSK
% curves are dashed same-colour and explained in title/caption).
nf_add_legend(fig7b, methods, styles6, colors6, ...
    'FontSize', 8, 'filled_idx', compressed_idx);
% Height 8.5 cm (+1.5 cm vs original 7.0)
nf_export_fig(gcf, 'fig7b_source_robustness', 'double', 'Height', 9.5);
fprintf('Fig.7b -> fig7b_source_robustness.pdf  |  CSV -> %s\n', CSV);
end

% ====================================================================
function res = mc_trial(Xm, Hm, tht, rt, Wm, Ym, Rm, P)
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
