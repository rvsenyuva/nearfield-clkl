function lh = nf_add_legend(fig, methods, styles, colors, varargin)
%NF_ADD_LEGEND  Place a single unified horizontal legend below all subplots.
%
%  lh = nf_add_legend(fig, methods, styles, colors)
%  lh = nf_add_legend(fig, methods, styles, colors, Name, Value, ...)
%
%  Compresses all data axes into a reserved interior region of the figure,
%  leaving BOTTOM_RESERVE at the bottom for the legend strip and TOP_RESERVE
%  at the top for the sgtitle, then places a single horizontal legend
%  centred at the bottom.
%
%  OPTIONAL NAME-VALUE PARAMETERS
%  -------------------------------
%  'FontSize'      : legend font size in pt        (default 8)
%  'filled_idx'    : indices of compressed-domain (filled-marker) methods
%  'NumColumns'    : number of legend columns       (default = numel(methods),
%                    i.e., single row; set to 4 for 8-entry 2-row legend)
%  'TopReserve'    : fraction of figure height kept free above subplots
%                    for sgtitle clearance           (default 0.05)
%
%  See also: nf_export_fig, plot_fig

p = inputParser();
p.addParameter('FontSize',    8,   @isnumeric);
p.addParameter('filled_idx',  [],  @isnumeric);
p.addParameter('NumColumns',  [],  @isnumeric);   % [] → all in one row
p.addParameter('TopReserve',  0.05, @isnumeric);  % space above subplots
p.parse(varargin{:});
opt = p.Results;

n_meth = numel(methods);
if isempty(opt.NumColumns)
    num_cols = n_meth;   % single-row default
else
    num_cols = opt.NumColumns;
end

% ── Step 1: compress data axes into the interior reserved region ──────────
% Bottom reserve: legend strip + x-axis tick label overhang.
% Top reserve:    sgtitle clearance (larger when subplot titles are multi-line).
BOTTOM_RESERVE = 0.15;
TOP_RESERVE    = opt.TopReserve;

drawnow;

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
    % Available interior height: from BOTTOM_RESERVE to (1 - TOP_RESERVE)
    interior_h = 1 - BOTTOM_RESERVE - TOP_RESERVE;
    new_bottom = pos(2) * interior_h + BOTTOM_RESERVE;
    new_height = pos(4) * interior_h;
    ax.Position = [pos(1), new_bottom, pos(3), new_height];
end

% ── Step 2: ghost axes at y=0 hosting dummy lines for the legend ─────────
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
        mfc = col;
    else
        mfc = 'none';
    end
    dum(mi) = plot(ax_leg, NaN, NaN, sty, ...
        'Color',           col, ...
        'MarkerFaceColor', mfc, ...
        'LineWidth',       1.5, ...
        'MarkerSize',      6, ...
        'DisplayName',     methods{mi});
end

% ── Step 3: create and position the legend ────────────────────────────────
lh = legend(ax_leg, dum, methods, ...
    'Orientation', 'horizontal', ...
    'NumColumns',  num_cols, ...
    'FontSize',    opt.FontSize, ...
    'Box',         'on');

drawnow;

lh.Units = 'normalized';
leg_w    = lh.Position(3);
leg_h    = lh.Position(4);
lh.Position = [(1 - leg_w) / 2,  0.01,  leg_w,  leg_h];
end
