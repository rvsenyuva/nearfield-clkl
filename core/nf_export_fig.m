function nf_export_fig(fig_handle, filename, layout, varargin)
%NF_EXPORT_FIG  Export a MATLAB figure to publication-quality PDF.
%
%  nf_export_fig(fig_handle, filename, layout)
%  nf_export_fig(fig_handle, filename, layout, Name, Value, ...)
%
%  Prepares a figure for inclusion in an IEEE two-column journal paper
%  (IEEEtran, 10pt body text) and exports it as a tight-cropped PDF.
%  All run_figN scripts call this function instead of plain saveas().
%
%  INPUTS
%  ------
%  fig_handle : figure handle (use gcf if called immediately after plotting)
%  filename   : output filename WITHOUT extension, e.g. 'fig2_nmse_snr'
%               The function appends '.pdf' automatically.
%               Output is written to the current working directory.
%  layout     : string, one of:
%                 'single'  -- one-column wide  (8.8 cm = 3.46 in)
%                 'double'  -- two-column wide  (18.2 cm = 7.17 in)  [default]
%                 'square'  -- equal width/height, one-column
%
%  OPTIONAL NAME-VALUE PAIRS
%  -------------------------
%  'Height'         : figure height in cm  (overrides layout default)
%  'FontSize'       : base axis font size in pt   (default 9)
%  'LegendFontSize' : legend font size in pt      (default 8)
%  'TitleFontSize'  : title font size in pt       (default 9)
%  'LabelFontSize'  : axis label font size in pt  (default 9)
%  'LineWidth'      : default line width in pt    (default 1.5)
%  'MarkerSize'     : default marker size in pt   (default 6)
%  'GridAlpha'      : grid transparency           (default 0.35)
%  'Verbose'        : print confirmation message  (default true)
%
%  OUTPUT
%  ------
%  Saves <filename>.pdf to the current working directory.
%  Returns nothing; the figure is left open for inspection.
%
%  USAGE EXAMPLES
%  --------------
%  % After running plot_fig() which calls saveas() already,
%  % replace the saveas() call with nf_export_fig():
%
%    plot_fig(SNR_vec, NMSE_db, NMSE_std, methods, ...
%        'Fig.2: NMSE vs SNR', 'SNR (dB)', 'Fig2');
%    nf_export_fig(gcf, 'fig2_nmse_snr', 'double');
%
%  % For a square convergence diagnostic figure (Fig. 9):
%    nf_export_fig(gcf, 'fig9_convergence', 'double', 'Height', 12);
%
%  % For Fig. 5 (near-to-far, wider aspect):
%    nf_export_fig(gcf, 'fig5_nearfar', 'double', 'Height', 9);
%
%  IEEEtran LAYOUT REFERENCE
%  --------------------------
%  Page width:        17.4 cm  (two 8.5 cm columns + 0.4 cm gutter)
%  Single column:      8.8 cm
%  Double column:     18.2 cm  (slight bleed into margins; safe for \figure*)
%  Body font size:    10 pt
%  Recommended axis font: 9 pt  (legible at final print scale)
%  Minimum axis font: 8 pt  (avoid going below this)
%
%  PDF vs EPS NOTE
%  ---------------
%  This function exports PDF because pdflatex handles PDF natively with
%  no conversion step.  EPS requires epstopdf or latex+dvipdf toolchain.
%  All run_figN scripts should use PDF export exclusively.
%
%  REQUIREMENTS
%  ------------
%  MATLAB R2020b or later.  No toolboxes required.
%  The export uses print() with -dpdf and -painters renderer for
%  vector output (no rasterisation of lines/text).
%
%  AXES SKIP TAGS
%  --------------
%  Axes whose Tag property matches any of the following are skipped
%  during font/line-width formatting (they are still exported):
%    'Colorbar'       -- colour-bar axes added by colorbar()
%    'legend'         -- MATLAB's internal legend axes
%    'legend_dummy'   -- ghost axes created by nf_add_legend()
%
%  See also: plot_fig, nf_add_legend, nf_params, run_fig2_nmse_snr

% ── Parse inputs ──────────────────────────────────────────────────────────
if nargin < 2 || isempty(filename)
    error('nf_export_fig: filename is required.');
end
if nargin < 3 || isempty(layout)
    layout = 'double';
end

% Strip extension if user accidentally included it
filename = regexprep(filename, '\.pdf$', '', 'ignorecase');

p = inputParser();
p.addParameter('Height',         [],    @isnumeric);
p.addParameter('FontSize',       9,     @isnumeric);
p.addParameter('LegendFontSize', 8,     @isnumeric);
p.addParameter('TitleFontSize',  9,     @isnumeric);
p.addParameter('LabelFontSize',  9,     @isnumeric);
p.addParameter('LineWidth',      1.5,   @isnumeric);
p.addParameter('MarkerSize',     6,     @isnumeric);
p.addParameter('GridAlpha',      0.35,  @isnumeric);
p.addParameter('Verbose',        true,  @islogical);
p.parse(varargin{:});
opt = p.Results;

% ── Column widths (cm) ────────────────────────────────────────────────────
switch lower(layout)
    case 'single'
        w_cm = 8.8;
        h_cm = 7.0;   % ~0.8 aspect ratio; readable in single column
    case 'double'
        w_cm = 18.2;
        h_cm = 7.5;   % matches Table* environment height in IEEEtran
    case 'square'
        w_cm = 8.8;
        h_cm = 8.8;
    otherwise
        warning('nf_export_fig: unknown layout "%s". Using ''double''.', layout);
        w_cm = 18.2;
        h_cm = 7.5;
end

% Override height if specified
if ~isempty(opt.Height)
    h_cm = opt.Height;
end

% Tags to skip during formatting (ghost/overlay axes that should not be
% touched by the font/linewidth pass below).
SKIP_TAGS = {'Colorbar', 'legend', 'legend_dummy'};

% ── Apply IEEE-standard formatting to all axes in the figure ─────────────
all_axes = findall(fig_handle, 'Type', 'axes');

for ax = all_axes'
    % Skip colour-bar axes, MATLAB's internal legend axes, and ghost axes
    % created by nf_add_legend().
    if isprop(ax, 'Tag') && any(strcmp(ax.Tag, SKIP_TAGS))
        continue
    end

    % Font sizes
    ax.FontSize  = opt.FontSize;
    ax.LabelFontSizeMultiplier = opt.LabelFontSize / opt.FontSize;
    ax.TitleFontSizeMultiplier = opt.TitleFontSize  / opt.FontSize;

    % Grid appearance
    ax.GridAlpha      = opt.GridAlpha;
    ax.MinorGridAlpha = opt.GridAlpha * 0.5;
    ax.Box            = 'off';   % cleaner look than box on in journal figs
    ax.TickDir        = 'out';   % standard IEEE journal tick direction

    % Apply line width and marker size to all Line objects in this axis
    lines_in_ax = findall(ax, 'Type', 'Line');
    for ln = lines_in_ax'
        if ln.LineWidth < opt.LineWidth
            ln.LineWidth = opt.LineWidth;
        end
        if ln.MarkerSize < opt.MarkerSize
            ln.MarkerSize = opt.MarkerSize;
        end
    end

    % ErrorBar objects (used in plot_fig)
    errorbars = findall(ax, 'Type', 'ErrorBar');
    for eb = errorbars'
        if eb.LineWidth < opt.LineWidth
            eb.LineWidth = opt.LineWidth;
        end
        if eb.MarkerSize < opt.MarkerSize
            eb.MarkerSize = opt.MarkerSize;
        end
    end
end

% Legend font size (applies to real legends, not the ghost axes)
all_legends = findall(fig_handle, 'Type', 'Legend');
for lg = all_legends'
    lg.FontSize = opt.LegendFontSize;
end

% ── Suppress axes toolbars (prevents MATLAB toolbar icons appearing in PDF) ─
all_ax_for_toolbar = findall(fig_handle, 'Type', 'axes');
for ax = all_ax_for_toolbar'
    try; ax.Toolbar.Visible = 'off'; catch; end
end

% ── Set figure paper size exactly to target dimensions ───────────────────
% Convert cm to inches for MATLAB paper units
w_in = w_cm / 2.54;
h_in = h_cm / 2.54;

set(fig_handle, ...
    'Units',           'inches', ...
    'Position',        [1 1 w_in h_in], ...
    'PaperUnits',      'inches', ...
    'PaperSize',       [w_in h_in], ...
    'PaperPosition',   [0 0 w_in h_in], ...
    'PaperPositionMode','manual', ...
    'Color',           'white');

% ── Export to PDF using vector renderer ──────────────────────────────────
outfile = [filename '.pdf'];
print(fig_handle, outfile, '-dpdf', '-painters', '-r0');

% ── Confirmation ─────────────────────────────────────────────────────────
if opt.Verbose
    fprintf('nf_export_fig: saved %s  [%.1f x %.1f cm, %s layout]\n', ...
        outfile, w_cm, h_cm, layout);
end
end
