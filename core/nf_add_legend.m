function lh = nf_add_legend(fig, methods, styles, colors, varargin)
%NF_ADD_LEGEND  Place a single unified horizontal legend below all subplots.
%
%  lh = nf_add_legend(fig, methods, styles, colors)
%  lh = nf_add_legend(fig, methods, styles, colors, 'FontSize', 8)
%
%  Compresses all data axes upward by BOTTOM_RESERVE (15% of figure height)
%  to create unobstructed space below the x-axis tick labels, then places
%  a single horizontal legend centred in that reserved strip.
%
%  The ghost axes used for the legend is tagged 'legend_dummy' so that
%  nf_export_fig skips it during font/linewidth formatting.
%
%  INPUTS
%  ------
%  fig     : figure handle
%  methods : cell array of method name strings
%  styles  : cell array of plot style strings
%  colors  : cell array of [R G B] colour vectors
%
%  OPTIONAL NAME-VALUE
%  -------------------
%  'FontSize'   : legend font size in pt  (default 8)
%  'filled_idx' : indices of compressed-domain (filled-marker) methods
%
%  See also: nf_export_fig, plot_fig

p = inputParser();
p.addParameter('FontSize',   8,  @isnumeric);
p.addParameter('filled_idx', [], @isnumeric);
p.addParameter('NumColumns', 0,  @isnumeric);  % 0=auto (all in one row)
p.parse(varargin{:});
opt = p.Results;

n_meth = numel(methods);

% ── Step 1: compress data axes upward to reserve bottom strip ────────────
% Maps each axis from figure-normalised [0,1] space into the top
% (1-BOTTOM_RESERVE) fraction, leaving BOTTOM_RESERVE free at the bottom
% for the legend strip and the downward-extending x-axis tick labels.
% Axes with height < 0.05 are skipped (sgtitle / annotation axes).
BOTTOM_RESERVE = 0.15;

drawnow;   % flush renderer so Position values are current

all_ax = findall(fig, 'Type', 'axes');
for ax = all_ax'
    if isprop(ax, 'Tag') && ...
       (strcmp(ax.Tag, 'Colorbar') || strcmp(ax.Tag, 'legend') || ...
        strcmp(ax.Tag, 'legend_dummy'))
        continue
    end
    pos = ax.Position;
    if pos(4) < 0.05
        continue   % skip sgtitle / annotation axes
    end
    % Remap bottom from [0,1] into [BOTTOM_RESERVE, 1]
    new_bottom = pos(2) * (1 - BOTTOM_RESERVE) + BOTTOM_RESERVE;
    new_height = pos(4) * (1 - BOTTOM_RESERVE);
    ax.Position = [pos(1), new_bottom, pos(3), new_height];
end

% ── Step 2: ghost axes at y=0 (invisible; hosts dummy lines for legend) ──
ax_leg = axes(fig, ...
    'Position',          [0, 0, 1, 0.001], ...
    'Visible',           'off', ...
    'Tag',               'legend_dummy', ...
    'HandleVisibility',  'off');
hold(ax_leg, 'on');

dum = gobjects(n_meth, 1);
for mi = 1:n_meth
    col = colors{mi};
    sty = styles{mi};
    if ~isempty(opt.filled_idx) && ismember(mi, opt.filled_idx)
        mfc = col;    % filled = compressed domain
    else
        mfc = 'none'; % open   = full-array
    end
    dum(mi) = plot(ax_leg, NaN, NaN, sty, ...
        'Color',           col, ...
        'MarkerFaceColor', mfc, ...
        'LineWidth',       1.5, ...
        'MarkerSize',      6, ...
        'DisplayName',     methods{mi});
end

% ── Step 3: create and position the legend ────────────────────────────────
n_cols_leg = opt.NumColumns;
if n_cols_leg <= 0
    n_cols_leg = numel(methods);  % all in one row (default)
end
lh = legend(ax_leg, dum, methods, ...
    'Orientation', 'horizontal', ...
    'NumColumns',  n_cols_leg, ...
    'FontSize',    opt.FontSize, ...
    'Box',         'on');

drawnow;   % flush so lh.Position reflects actual rendered legend size

lh.Units = 'normalized';
leg_w    = lh.Position(3);
leg_h    = lh.Position(4);

% Centre horizontally; clamp so legend never clips either edge.
x_left = max(0.01, (1 - leg_w) / 2);
if x_left + leg_w > 0.99
    x_left = max(0.01, 0.99 - leg_w);
end
lh.Position = [x_left,  0.01,  leg_w,  leg_h];
end
