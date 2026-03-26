function save_results(tag, NMSE_db, NMSE_std, x_vec, methods, varargin)
%SAVE_RESULTS  Save simulation results to a .mat file for reproducibility.
%
%  save_results(tag, NMSE_db, NMSE_std, x_vec, methods)
%  save_results(tag, NMSE_db, NMSE_std, x_vec, methods, 'extra_key', extra_val, ...)
%
%  Creates <tag>.mat in the current directory with all inputs plus a
%  timestamp and a struct of optional extra arrays (e.g., RMSE_r, RT_mean).

S.tag       = tag;
S.NMSE_db   = NMSE_db;
S.NMSE_std  = NMSE_std;
S.x_vec     = x_vec;
S.methods   = methods;
S.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% Parse optional key-value pairs
for k = 1:2:numel(varargin)
    key = varargin{k};
    val = varargin{k+1};
    S.(key) = val;
end

fname = [tag, '.mat'];
save(fname, '-struct', 'S');
fprintf('Results saved to %s\n', fname);
end
