function run_all(mode)
%RUN_ALL  Master runner -- executes all ten figures (Steps 1-7 of TWC plan).
%
%  USAGE
%  -----
%    run_all                    % full run, serial       (N_MC=400, ~2-4 h)
%    run_all('fast')            % preview, serial        (N_MC=20,  ~5 min)
%    run_all('parfor')          % full run, parallel     (N_MC=400, ~10-20 min)
%    run_all('fast parfor')     % preview, parallel      (N_MC=20,  ~1 min)
%    run_all('parfor fast')     % same -- keyword order does not matter
%
%  PARALLELISM NOTES
%  -----------------
%  The 'parfor' flag distributes the Monte Carlo trials in each sweep
%  point across parallel workers.  Requirements:
%    _ MATLAB Parallel Computing Toolbox must be installed.
%    _ A parpool is opened automatically on first parfor call (~10-30 s).
%    _ To control worker count: parpool('local', N) before calling run_all.
%    _ Fig.6 (runtime) always runs serially -- parfor distorts wall-clock
%      timings that are the purpose of that figure.
%  If the Parallel Computing Toolbox is absent, parfor silently falls back
%  to a serial for-loop so the code is always safe to run.
%
%  DESIGN
%  ------
%  Each run_fig*.m function is a proper MATLAB function (not a script).
%  This gives it an isolated workspace so internal 'clear' calls cannot
%  corrupt run_all's local state (n_ok, t_total, etc.).

if nargin < 1; mode = 'full'; end
mode_lower = lower(mode);
fast    = contains(mode_lower, 'fast');
use_par = contains(mode_lower, 'parfor');

if fast
    fprintf('\n*** FAST MODE: N_MC = 20 in every figure ***\n\n');
end
if use_par
    fprintf('*** PARALLEL MODE: parfor enabled across MC trials ***\n');
    if isempty(ver('parallel'))
        warning('run_all:noPCT', ...
            'Parallel Computing Toolbox not found -- running serially.');
        use_par = false;
    else
        % Open a parpool if one is not already running
        if isempty(gcp('nocreate'))
            fprintf('  Opening parpool ...\n');
            parpool('local');
        end
        p = gcp('nocreate');
        fprintf('  Workers: %d\n\n', p.NumWorkers);
    end
end

% ---- Reproducibility -----------------------------------------------
% Reset global RNG here so run_all always produces the same numbers.
% Each run_fig*.m also resets the same seed for standalone use.
rng(42, 'twister');  % Fixed seed: ensures bit-exact reproducibility across runs
if exist('nf_simulation_results.csv','file')
    delete('nf_simulation_results.csv');
    fprintf('  Cleared nf_simulation_results.csv for fresh run.\n');
end

print_table1();

t_total = tic;
n_ok    = 0;

n_ok = run_one('Fig.2  NMSE vs SNR',   @() run_fig2_nmse_snr(fast, use_par),   n_ok);
n_ok = run_one('Fig.3  NMSE vs N_RF',  @() run_fig3_nmse_nrf(fast, use_par),   n_ok);
n_ok = run_one('Fig.4  NMSE vs N',     @() run_fig4_nmse_n(fast, use_par),     n_ok);
n_ok = run_one('Fig.5  Near to Far',   @() run_fig5_nearfar(fast, use_par),    n_ok);
n_ok = run_one('Fig.6  Runtime vs M',  @() run_fig6_runtime(fast),             n_ok);  % always serial
n_ok = run_one('Fig.7  Robustness',    @() run_fig7_robustness(fast, use_par), n_ok);
n_ok = run_one('Fig.8  NMSE vs d',     @() run_fig8_vard(fast, use_par),       n_ok);
n_ok = run_one('Fig.9  Convergence',   @() run_fig9_convergence(fast),         n_ok);
n_ok = run_one('Fig.7b Source robust', @() run_fig7b_source_robustness(fast, use_par), n_ok);

fprintf('\n%s\n  Done: %d/10 figures  |  %.1f min total\n%s\n\n', ...
    repmat('=',1,55), n_ok, toc(t_total)/60, repmat('=',1,55));
end


% -------------------------------------------------------------------
function n_ok = run_one(label, fn, n_ok)
fprintf('\n%s\n  Running: %s\n%s\n', repmat('-',1,55), label, repmat('-',1,55));
t0 = tic;
try
    fn();
    fprintf('  --> OK  (%.1f s)\n', toc(t0));
    n_ok = n_ok + 1;
catch ME
    fprintf('  --> FAILED: %s\n  at %s (line %d)\n', ...
        ME.message, ME.stack(1).name, ME.stack(1).line);
end
end


% -------------------------------------------------------------------
function print_table1()
fprintf('\n%s\n  TABLE I: Complexity Comparison\n%s\n', ...
    repmat('-',1,90), repmat('-',1,90));
fprintf('%-26s | %-22s | %-28s | %s\n', ...
    'Method','Dict. size','Dominant cost','Notes');
fprintf('%s\n', repmat('-',1,90));
rows = { ...
  'Polar-grid SOMP',     'Q_th x Q_r atoms',   'O(Q_th*Q_r*N_RF^2)/iter', 'Coherence grows with Q_r'; ...
  'DL-OMP (Zhang)',      'NxN dict.',           'O(3K_iter*P*N)',           'Zhang et al. 2024, full-array'; ...
  'DFrFT-NOMP',          '(gridless)',          'O(PM log M_)+Newton',     'Yang et al. 2024, full-array'; ...
  'MUSIC+Triangulation', '(no dict.)',          'O(Q*M_s^3)',              'Low-complexity ref.'; ...
  'Proposed CL-KL',      'Q_th (angle only)',   'O(N_RF^3+Q_th*N_RF^2)',   'Multi-start, ring-indexed init' ...
};
for k = 1:size(rows,1)
    fprintf('%-26s | %-22s | %-28s | %s\n', rows{k,:});
end
fprintf('%s\n\n', repmat('-',1,90));
end
