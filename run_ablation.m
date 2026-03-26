function run_ablation(use_par)
%RUN_ABLATION  Ablation study for CL-KL (P2 response to reviewer Q1).
%
%  Quantifies the individual contribution of each CL-KL design choice:
%    (a) Frozen N0 estimate (vs. iterative residual re-estimation)
%    (b) Multi-start warm-start (vs. single far-range start only)
%    (c) Post-loop 4-pass MF scan (vs. warm-start u as final estimate)
%
%  4 configurations x 3 SNR points x N_MC=50 trials.
%  Reports NMSE (dB), RMSE_theta (deg), RMSE_r (m) for each cell.
%
%  Usage:
%    run_ablation()          -- serial (default)
%    run_ablation(true)      -- parfor across MC trials
%
%  Output files:
%    ablation_results.mat         -- full struct for post-processing
%    nf_ablation_results.csv      -- 9-column CSV (one row per config x SNR)
%
%  NOTE: All ablation flags in nf_clkl.m default to false; this script
%  is the ONLY place they are ever set to true.  Production run scripts
%  are unaffected.

if nargin < 1; use_par = false; end

% ---- Setup ----------------------------------------------------------
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);

% ---- Parameters (standard Fig.2 configuration) ----------------------
P      = nf_params();
P.M    = 64;
P      = nf_update_derived_pub(P);
P.N_RF = 8;
P.N    = 64;
P.d    = 3;

N_MC   = 50;       % sufficient for ablation; not a publication figure
P.N_MC = N_MC;

SNR_dB_vec = [-5, 5, 15];   % low / mid / high SNR
n_snr      = numel(SNR_dB_vec);

% ---- Ablation configurations ----------------------------------------
%  Config 1: Full CL-KL  -- reference, all flags false
%  Config 2: w/o frozen N0  -- re-estimate N0 from residual every 10 iters
%  Config 3: w/o multi-start  -- only far-range init (Start 3), skip ring/near
%  Config 4: w/o post-loop scan  -- skip Phase 2, use warm-start u as final
cfg_names = {'CL-KL (full)', 'w/o frozen N0', 'w/o multi-start', 'w/o post-loop scan'};
cfg_N0    = [false, true,  false, false];
cfg_1st   = [false, false, true,  false];
cfg_ns    = [false, false, false, true ];
n_cfg     = numel(cfg_names);

% ---- Result storage -------------------------------------------------
NMSE_dB  = nan(n_cfg, n_snr);
RMSE_th  = nan(n_cfg, n_snr);   % degrees
RMSE_r   = nan(n_cfg, n_snr);   % metres
fail_pct = nan(n_cfg, n_snr);   % percent

% ---- Monte Carlo loop -----------------------------------------------
fprintf('=== Ablation Study: %d configs x %d SNR points x N_MC=%d ===\n\n', ...
    n_cfg, n_snr, N_MC);

for si = 1:n_snr
    SNR = SNR_dB_vec(si);
    fprintf('SNR = %+3d dB  -- pre-generating %d trials ... ', SNR, N_MC);

    % Pre-generate all N_MC trials once (same data for all 4 configs)
    [X_all, H_all, th_all, r_all, W_all, Y_all, R_all] = pregen(P, SNR, true);
    fprintf('done.\n');

    for ci = 1:n_cfg
        % Build per-config P struct (only ablation flags differ)
        P_run = P;
        P_run.ablation_update_N0    = cfg_N0(ci);
        P_run.ablation_single_start = cfg_1st(ci);
        P_run.ablation_skip_scan    = cfg_ns(ci);

        nmse_mc = nan(1, N_MC);
        rth_mc  = nan(1, N_MC);
        rr_mc   = nan(1, N_MC);
        fl_mc   = nan(1, N_MC);

        if use_par
            parfor mc = 1:N_MC
                try
                    [th_h, r_h] = nf_clkl(R_all{mc}, W_all{mc}, P_run);
                    [nm, rth, rr, fl] = nf_metrics(th_h, r_h, Y_all{mc}, ...
                        W_all{mc}, H_all{mc}, th_all{mc}, r_all{mc}, P_run);
                    nmse_mc(mc) = nm;
                    rth_mc(mc)  = rth;
                    rr_mc(mc)   = rr;
                    fl_mc(mc)   = fl;
                catch
                    nmse_mc(mc) = 1.0; rth_mc(mc) = pi;
                    rr_mc(mc)   = 1e6; fl_mc(mc)  = 1.0;
                end
            end
        else
            for mc = 1:N_MC
                try
                    [th_h, r_h] = nf_clkl(R_all{mc}, W_all{mc}, P_run);
                    [nmse_mc(mc), rth_mc(mc), rr_mc(mc), fl_mc(mc)] = ...
                        nf_metrics(th_h, r_h, Y_all{mc}, W_all{mc}, ...
                                   H_all{mc}, th_all{mc}, r_all{mc}, P_run);
                catch
                    nmse_mc(mc) = 1.0; rth_mc(mc) = pi;
                    rr_mc(mc)   = 1e6; fl_mc(mc)  = 1.0;
                end
            end
        end

        % Aggregate: nanmean excludes any remaining NaN trials
        NMSE_dB(ci, si)  = 10*log10(nanmean(nmse_mc));
        RMSE_th(ci, si)  = nanmean(rth_mc) * 180/pi;   % rad -> deg
        RMSE_r(ci, si)   = nanmean(rr_mc);
        fail_pct(ci, si) = 100 * nanmean(fl_mc);

        fprintf('  [%d/4] %-24s  NMSE=%6.2f dB  RMSE_th=%5.3f deg  RMSE_r=%5.1f m  fail=%4.1f%%\n', ...
            ci, cfg_names{ci}, NMSE_dB(ci,si), RMSE_th(ci,si), RMSE_r(ci,si), fail_pct(ci,si));
    end
    fprintf('\n');
end

% ---- Print summary table -------------------------------------------
print_table(cfg_names, SNR_dB_vec, NMSE_dB, RMSE_th, RMSE_r, fail_pct);

% ---- Save results --------------------------------------------------
save('ablation_results.mat', 'cfg_names', 'SNR_dB_vec', ...
     'NMSE_dB', 'RMSE_th', 'RMSE_r', 'fail_pct', 'N_MC', 'P');
fprintf('\nResults saved to ablation_results.mat\n');

% ---- Write nf_ablation_results.csv ---------------------------------
% Schema (9 columns):
%   config_name, SNR_dB, NMSE_dB_mean, RMSE_theta_deg_mean,
%   RMSE_r_m_mean, fail_rate_pct, N_MC, M, N_RF
save_ablation_csv('nf_ablation_results.csv', cfg_names, SNR_dB_vec, ...
    NMSE_dB, RMSE_th, RMSE_r, fail_pct, N_MC, P);
fprintf('Results saved to nf_ablation_results.csv\n');
end


% ====================================================================
function print_table(names, snr_vec, NMSE_dB, RMSE_th, RMSE_r, fail_pct)
%PRINT_TABLE  Console-formatted summary table for all three metrics.
n_cfg = numel(names);
n_snr = numel(snr_vec);

div = repmat('=', 1, 76);
sep = repmat('-', 1, 76);

fprintf('\n%s\n', div);
fprintf('  ABLATION RESULTS  (N_MC=50, M=64, N_RF=8, N=64, d=3)\n');
fprintf('%s\n\n', div);

% Build SNR header string
hdr = sprintf('%-26s', 'Configuration');
for si = 1:n_snr
    hdr = [hdr, sprintf(' %+3ddB ', snr_vec(si))]; %#ok<AGROW>
end

% ----- NMSE (dB) -----
fprintf('CHANNEL NMSE (dB):\n%s\n', sep);
fprintf('%s\n', hdr);
fprintf('%s\n', sep);
for ci = 1:n_cfg
    row = sprintf('%-26s', names{ci});
    for si = 1:n_snr
        row = [row, sprintf(' %5.2f  ', NMSE_dB(ci,si))]; %#ok<AGROW>
    end
    if ci == 1
        row = [row, ' <- reference']; %#ok<AGROW>
    else
        % Degradation relative to Config 1 at mid SNR (SNR index 2)
        delta = NMSE_dB(ci,2) - NMSE_dB(1,2);
        row = [row, sprintf(' (+%.2f dB vs full at SNR=%+ddB)', delta, snr_vec(2))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('\n');

% ----- RMSE_theta (deg) -----
fprintf('ANGLE RMSE (deg):\n%s\n', sep);
fprintf('%s\n', hdr);
fprintf('%s\n', sep);
for ci = 1:n_cfg
    row = sprintf('%-26s', names{ci});
    for si = 1:n_snr
        row = [row, sprintf(' %5.3f  ', RMSE_th(ci,si))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('\n');

% ----- RMSE_r (m) -----
fprintf('RANGE RMSE (m):\n%s\n', sep);
fprintf('%s\n', hdr);
fprintf('%s\n', sep);
for ci = 1:n_cfg
    row = sprintf('%-26s', names{ci});
    for si = 1:n_snr
        row = [row, sprintf(' %6.1f  ', RMSE_r(ci,si))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('\n');

% ----- Failure rate (%) -----
fprintf('FAILURE RATE (%%):\n%s\n', sep);
fprintf('%s\n', hdr);
fprintf('%s\n', sep);
for ci = 1:n_cfg
    row = sprintf('%-26s', names{ci});
    for si = 1:n_snr
        row = [row, sprintf('  %4.1f   ', fail_pct(ci,si))]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('%s\n\n', div);
end


% ====================================================================
function save_ablation_csv(csv_path, cfg_names, snr_vec, ...
    NMSE_dB, RMSE_th, RMSE_r, fail_pct, N_MC, P)
%SAVE_ABLATION_CSV  Write ablation results to a 9-column CSV.
%
%  One row per (config x SNR) combination = 4 configs x 3 SNR = 12 rows.
%  Appends timestamp so successive runs accumulate.

HEADER = 'timestamp,config_name,ablation_flags,SNR_dB,NMSE_dB_mean,RMSE_theta_deg_mean,RMSE_r_m_mean,fail_rate_pct,N_MC';

% Ablation flag descriptions per config
flags = {'none (full CL-KL)', 'update_N0', 'single_start', 'skip_scan'};

write_header = ~isfile(csv_path);
fid = fopen(csv_path, 'a');
if write_header
    fprintf(fid, '%s\n', HEADER);
end

ts = datestr(now, 'yyyy-mm-dd HH:MM:SS');
for ci = 1:numel(cfg_names)
    for si = 1:numel(snr_vec)
        fprintf(fid, '%s,%s,%s,%d,%.4f,%.4f,%.4f,%.2f,%d\n', ...
            ts, cfg_names{ci}, flags{ci}, snr_vec(si), ...
            NMSE_dB(ci,si), RMSE_th(ci,si), RMSE_r(ci,si), ...
            fail_pct(ci,si), N_MC);
    end
end
fclose(fid);
end
