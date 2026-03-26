function plot_fig(x_vec, NMSE_db, NMSE_std, methods, title_str, xlabel_str, fig_name)
%PLOT_FIG  Standard errorbar plot used by Figs 2, 3, 4, 7.
%
%  plot_fig(x_vec, NMSE_db, NMSE_std, methods, title_str, xlabel_str, fig_name)
%
%  CHANGE (legend fix): legend placed below axes ('southoutside', horizontal,
%  3 columns) so it never obstructs plotted curves.  Callers must increase
%  the nf_export_fig Height by ~2 cm to accommodate the extra legend strip.

n_meth = numel(methods);
% Method order: CL-KL, P-SOMP, DL-OMP, MUSIC+Tri, DFrFT-NOMP, BF-SOMP
styles = {'-o','-s','--^','-.d','-v','-.h'};
colors = {[0 0.45 0.74],[0.85 0.33 0.10],[0.47 0.67 0.19], ...
          [0.49 0.18 0.56],[0.93 0.69 0.13],[0.30 0.57 0.43]};

% Compressed-domain methods (filled markers): CL-KL, P-SOMP
% Full-array methods (open markers): DL-OMP, MUSIC+Tri, DFrFT-NOMP, BF-SOMP
% This visual distinction encodes the 8x data volume difference (N_RF^2 vs M*N).
compressed_idx = [1, 2];   % CL-KL, P-SOMP use N_RF x N_RF compressed R_hat

figure('Name', fig_name, 'Position', [100 100 580 440]);
for mi = 1:n_meth
    col = colors{mod(mi-1,numel(colors))+1};
    if ismember(mi, compressed_idx)
        mfc = col;     % filled marker = compressed domain
    else
        mfc = 'none';  % open marker   = full-array access
    end
    errorbar(x_vec, NMSE_db(mi,:), NMSE_std(mi,:)/2, ...
        styles{mod(mi-1,numel(styles))+1}, ...
        'Color', col, 'MarkerFaceColor', mfc, ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', methods{mi});
    hold on;
end
xlabel(xlabel_str, 'FontSize', 12);
ylabel('NMSE (dB)', 'FontSize', 12);
title(title_str, 'FontSize', 13);
grid on;
set(gca, 'FontSize', 11);

% Legend below axes -- no curve obstruction.
% NumColumns=3 keeps the strip compact (2 rows of 3 for 6 methods).
leg = legend('Orientation', 'horizontal', ...
             'NumColumns',  3, ...
             'Location',    'southoutside', ...
             'FontSize',    9);
leg.Title.String    = 'Filled: compressed | Open: full-array';
leg.Title.FontSize  = 8;
end
