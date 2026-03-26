%SETUP  Add the simulation directory AND all subfolders to MATLAB path.
%
%  Run this once per MATLAB session before any other script:
%    >> cd('/path/to/your/nearfield_sim_folder')
%    >> setup
%  Or run it from anywhere with the full path:
%    >> run('/path/to/nearfield_sim_folder/setup.m')
%
%  Uses genpath() so that subfolders (core/, estimators/, figures/,
%  metrics/) are all added automatically, regardless of whether the
%  .m files are flat in one folder or organised into subdirectories.

this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(this_dir));
fprintf('Path configured (with subfolders): %s\n', this_dir);
fprintf('Ready to run: test_modules, test_fixes, run_all\n');