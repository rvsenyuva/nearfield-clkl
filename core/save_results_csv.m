function save_results_csv(csv_path, fig_id, fig_name, experiment_type, ...
        fixed_M, fixed_N_RF, fixed_N, fixed_d, fixed_SNR_dB, ...
        sweep_var, sweep_val, truth_model, ...
        methods, NMSE_dB_mean, NMSE_dB_std, ...
        RMSE_theta_deg, RMSE_r_m, fail_rate_pct, ...
        runtime_s, notes_cell, N_MC, r_RD_m, ...
        clkl_iters_mean, clkl_conv_pct, clkl_N0_mean, ...
        crb_theta_deg, crb_r_m_val)
%SAVE_RESULTS_CSV  Append one row per method to the master results CSV.
%
%  Called after each sweep point -- partial runs are preserved.
%  Creates the file with a header on first call; appends thereafter.
%
%  NEW COLUMNS vs previous version:
%    N_MC            -- trials used (20 fast / 200 full)
%    r_RD_m          -- Rayleigh distance [m]
%    clkl_avg_iters  -- mean CL-KL iterations (NaN for other methods)
%    clkl_conv_pct   -- % trials CL-KL converged (NaN for other methods)
%    clkl_N0_mean    -- mean final CL-KL N0 estimate (NaN for other methods)
%
%  The three clkl_* columns are filled only for the CL-KL row;
%  all other method rows get blank cells.  This makes it easy to
%  filter in the CSV: WHERE method = 'CL-KL' AND clkl_conv_pct < 50
%  immediately flags non-convergence problems.
%
%  CSV COLUMNS (25 total)
%  timestamp, figure_id, figure_name, experiment_type,
%  fixed_M, fixed_N_RF, fixed_N, fixed_d, fixed_SNR_dB,
%  sweep_variable, sweep_value, truth_model,
%  N_MC, r_RD_m,
%  method,
%  NMSE_dB_mean, NMSE_dB_std,
%  RMSE_theta_deg_mean, RMSE_r_m_mean, fail_rate_pct,
%  runtime_s,
%  clkl_avg_iters, clkl_conv_pct, clkl_N0_mean,
%  notes

% CSV COLUMNS (28 total -- 3 new vs v1: input_data, CRB_theta_deg, CRB_r_m)
HEADER = ['timestamp,figure_id,figure_name,experiment_type,' ...
          'fixed_M,fixed_N_RF,fixed_N,fixed_d,fixed_SNR_dB,' ...
          'sweep_variable,sweep_value,truth_model,' ...
          'N_MC,r_RD_m,' ...
          'method,input_data,' ...
          'NMSE_dB_mean,NMSE_dB_std,' ...
          'RMSE_theta_deg_mean,RMSE_r_m_mean,fail_rate_pct,' ...
          'runtime_s,' ...
          'clkl_avg_iters,clkl_conv_pct,clkl_N0_mean,' ...
          'CRB_theta_deg,CRB_r_m,' ...
          'notes'];

n_meth = numel(methods);

% ---- Handle optional trailing arguments ----------------------------
if nargin < 21 || isempty(N_MC);           N_MC           = NaN; end
if nargin < 22 || isempty(r_RD_m);         r_RD_m         = NaN; end
if nargin < 23 || isempty(clkl_iters_mean);clkl_iters_mean= NaN; end
if nargin < 24 || isempty(clkl_conv_pct);  clkl_conv_pct  = NaN; end
if nargin < 25 || isempty(clkl_N0_mean);   clkl_N0_mean   = NaN; end

if nargin < 26; crb_theta_deg = NaN; end
if nargin < 27; crb_r_m_val   = NaN; end

% ---- Pad per-method vectors ----------------------------------------
NMSE_dB_std    = pad_vec(NMSE_dB_std,    n_meth);
RMSE_theta_deg = pad_vec(RMSE_theta_deg, n_meth);
RMSE_r_m       = pad_vec(RMSE_r_m,       n_meth);
fail_rate_pct  = pad_vec(fail_rate_pct,  n_meth);
runtime_s      = pad_vec(runtime_s,      n_meth);
if isempty(notes_cell) || ~iscell(notes_cell) || numel(notes_cell) ~= n_meth
    notes_cell = repmat({''}, n_meth, 1);
end

% ---- Open file (append; write header if new) -----------------------
need_header = ~isfile(csv_path);
fid = fopen(csv_path, 'a');
if fid == -1
    warning('save_results_csv: cannot open %s for writing.', csv_path);
    return
end
if need_header
    fprintf(fid, '%s\n', HEADER);
end

% ---- input_data: classify each method as compressed or full-array --
% Compressed methods observe R_hat (N_RF x N_RF); others observe X (M x N)
input_data_map = containers.Map(...
    {'CL-KL','P-SOMP','DL-OMP','MUSIC+Tri','DFrFT-NOMP','BF-SOMP','BF-SOMP (adpt)'}, ...
    {'R_hat','R_hat','X','X','X','X','X'});

% ---- One row per method --------------------------------------------
ts = datestr(now, 'yyyy-mm-dd HH:MM:SS');
for mi = 1:n_meth
    % CL-KL diagnostics only apply to method index 1 (CL-KL row)
    % All other methods get blank cells for those columns
    if mi == 1
        c_iters = fmt(clkl_iters_mean);
        c_conv  = fmt(clkl_conv_pct);
        c_N0    = fmt(clkl_N0_mean);
    else
        c_iters = '';
        c_conv  = '';
        c_N0    = '';
    end

    % Resolve input_data string for this method
    if isKey(input_data_map, methods{mi})
        inp_data = input_data_map(methods{mi});
    else
        inp_data = 'X';  % default: full-array
    end
    % CRB and CL-KL diag columns only in CL-KL row (mi==1)
    if mi == 1
        c_crb_th = fmt(crb_theta_deg);
        c_crb_r  = fmt(crb_r_m_val);
    else
        c_crb_th = '';
        c_crb_r  = '';
    end
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
        ts,                        ...  % 1  timestamp
        fig_id,                    ...  % 2  figure_id
        fig_name,                  ...  % 3  figure_name
        experiment_type,           ...  % 4  experiment_type
        num2str(fixed_M),          ...  % 5  fixed_M
        fmt(fixed_N_RF),           ...  % 6  fixed_N_RF
        fmt(fixed_N),              ...  % 7  fixed_N
        fixed_d,                   ...  % 8  fixed_d
        fmt(fixed_SNR_dB),         ...  % 9  fixed_SNR_dB
        sweep_var,                 ...  % 10 sweep_variable
        fmt(sweep_val),            ...  % 11 sweep_value
        truth_model,               ...  % 12 truth_model
        fmt(N_MC),                 ...  % 13 N_MC
        fmt(r_RD_m),               ...  % 14 r_RD_m
        methods{mi},               ...  % 15 method
        inp_data,                  ...  % 16 input_data  (NEW)
        fmt(NMSE_dB_mean(mi)),     ...  % 17 NMSE_dB_mean
        fmt(NMSE_dB_std(mi)),      ...  % 18 NMSE_dB_std
        fmt(RMSE_theta_deg(mi)),   ...  % 19 RMSE_theta_deg_mean
        fmt(RMSE_r_m(mi)),         ...  % 20 RMSE_r_m_mean
        fmt(fail_rate_pct(mi)),    ...  % 21 fail_rate_pct
        fmt(runtime_s(mi)),        ...  % 22 runtime_s
        c_iters,                   ...  % 23 clkl_avg_iters
        c_conv,                    ...  % 24 clkl_conv_pct
        c_N0,                      ...  % 25 clkl_N0_mean
        c_crb_th,                  ...  % 26 CRB_theta_deg (NEW)
        c_crb_r,                   ...  % 27 CRB_r_m       (NEW)
        notes_cell{mi});               % 28 notes
end
fclose(fid);
end


function v = pad_vec(v, n)
if isempty(v) || numel(v) ~= n
    v = nan(n,1);
else
    v = v(:);
end
end

function s = fmt(x)
if isempty(x) || (isscalar(x) && isnan(x))
    s = '';
else
    s = sprintf('%.6g', x);
end
end
