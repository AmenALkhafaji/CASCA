%% ============================================================
% SCHUBERT PANEL DATA:
% C-ASCA LOADINGS + METHOD-WISE EFFECT HEATMAP
%
% Panel A = Factor 1 / method-related score plot
% Panel B = Factor 2 / group-related score plot
% Panel C = Combined model score plot
% Panel D = C-ASCA loadings: Disagreement + Consensus
% Panel E = Median(Case) - Median(Control), NOT oMEDA
%
% Final-panel fix:
%   - D and E have identical vertical heatmap size
%   - E stars are moved to upper part of each cell
%   - E values are slightly smaller to avoid star overlap
%   - Genus panel height increased slightly for clean 25-row display
%% ============================================================

clc;

%% ============================================================
% REQUIRED VARIABLES IN WORKSPACE:
%   X
%   Group
%   lab
%   ascao
%
% Optional:
%   outDir
%   exportResolution
%   analysisLevel
%% ============================================================

%% ============================================================
% OUTPUT SETTINGS
%% ============================================================

if ~exist('outDir', 'var') || isempty(outDir)
    outDir = fullfile(pwd, 'Figures_JPG');
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if ~exist('exportResolution', 'var') || isempty(exportResolution)
    exportResolution = 600;
end

if ~exist('analysisLevel', 'var') || isempty(analysisLevel)
    analysisLevel = "genus";
end

alpha = 0.05;
nDec = 2;
topK = 25;

taxaNames = string(lab(:));
nTaxa = numel(taxaNames);

%% ============================================================
% METHOD AND GROUP LABELS FOR SCORE PANELS
%% ============================================================

methodOrder = ["Raw Data","TSS","Rarefaction","CLR"];
groupOrder  = ["Control","Case"];

nMethods = numel(methodOrder);

if mod(size(X,1), nMethods) ~= 0
    error('Number of rows in X must be divisible by number of methods.');
end

N = size(X,1) / nMethods;

Mtr = strings(size(X,1),1);
Mtr(1:N)         = "Raw Data";
Mtr(N+1:2*N)     = "TSS";
Mtr(2*N+1:3*N)   = "Rarefaction";
Mtr(3*N+1:4*N)   = "CLR";

G = string(Group);
G(G == "-1") = "Control";
G(G == "1")  = "Case";

%% ============================================================
% SPLIT METHOD BLOCKS
%% ============================================================

idxRAW  = 1:N;
idxTSS  = N+1 : 2*N;
idxRARE = 2*N+1 : 3*N;
idxCLR  = 3*N+1 : 4*N;

X_RAW  = X(idxRAW,  :);
X_TSS  = X(idxTSS,  :);
X_RARE = X(idxRARE, :);
X_CLR  = X(idxCLR,  :);

g_RAW  = Group(idxRAW);
g_TSS  = Group(idxTSS);
g_RARE = Group(idxRARE);
g_CLR  = Group(idxCLR);

%% ============================================================
% BUILD C-ASCA CONSENSUS REPRESENTATION
%
% C-ASCA consensus representation = group/consensus + observation + depth
% If factors 3 and 4 are unavailable, fallback to factor 2 only.
%% ============================================================

xB = ascao.factors{2}.matrix;

if numel(ascao.factors) >= 4
    xC = ascao.factors{3}.matrix;
    xD = ascao.factors{4}.matrix;
    Xmb_keep = xB + xC + xD;
else
    warning('ascao.factors{3} or ascao.factors{4} not available. Using factor 2 only for C-ASCA effect heatmap.');
    Xmb_keep = xB;
end

if mod(size(Xmb_keep,1), N) ~= 0
    error('Rows of C-ASCA consensus matrix are not divisible by N.');
end

K = size(Xmb_keep,1) / N;
K = round(K);

X3 = reshape(Xmb_keep, N, K, []);
X_CASCA = squeeze(mean(X3, 2, 'omitnan'));

if size(X_CASCA,1) ~= N || size(X_CASCA,2) ~= nTaxa
    error('X_CASCA must be N x nTaxa.');
end

g_CASCA = g_RAW;

%% ============================================================
% FEATURE-WISE KRUSKAL-WALLIS + BH-FDR
%
% Effect direction:
%   positive = higher in Case
%   negative = higher in Control
%% ============================================================

[p_RAW,   q_RAW,   sig_RAW,   eff_RAW]   = local_kw_bh_feature_test(X_RAW,   g_RAW,   alpha);
[p_TSS,   q_TSS,   sig_TSS,   eff_TSS]   = local_kw_bh_feature_test(X_TSS,   g_TSS,   alpha);
[p_RARE,  q_RARE,  sig_RARE,  eff_RARE]  = local_kw_bh_feature_test(X_RARE,  g_RARE,  alpha);
[p_CLR,   q_CLR,   sig_CLR,   eff_CLR]   = local_kw_bh_feature_test(X_CLR,   g_CLR,   alpha);
[p_CASCA, q_CASCA, sig_CASCA, eff_CASCA] = local_kw_bh_feature_test(X_CASCA, g_CASCA, alpha);

Effect_all = [ ...
    eff_RAW(:), ...
    eff_TSS(:), ...
    eff_RARE(:), ...
    eff_CLR(:), ...
    eff_CASCA(:)];

Sig_all = [ ...
    sig_RAW(:), ...
    sig_TSS(:), ...
    sig_RARE(:), ...
    sig_CLR(:), ...
    sig_CASCA(:)];

methodNamesEffect = ["Raw Data","TSS","Rarefaction","CLR","C-ASCA"];

%% ============================================================
% C-ASCA LOADINGS
%
% Consensus direction is aligned with C-ASCA effect:
% positive consensus = higher in Case
% negative consensus = higher in Control
%% ============================================================

load_disagree_raw  = -ascao.factors{1}.loads(:,1);
load_consensus_raw = -ascao.factors{2}.loads(:,1);

if numel(load_disagree_raw) ~= nTaxa || numel(load_consensus_raw) ~= nTaxa
    error('ASCA loading length does not match number of taxa.');
end

refEffect = eff_CASCA(:);

rAlign = corr(load_consensus_raw(:), refEffect(:), 'Rows', 'complete');

if isnan(rAlign)
    warning('Consensus loading alignment correlation is NaN. Keeping original consensus direction.');
    load_consensus = load_consensus_raw;
elseif rAlign < 0
    load_consensus = -load_consensus_raw;
    fprintf('Consensus loading flipped to align with Median(Case) - Median(Control).\n');
else
    load_consensus = load_consensus_raw;
    fprintf('Consensus loading already aligned with Median(Case) - Median(Control).\n');
end

load_disagree = load_disagree_raw;

%% ============================================================
% SELECT TOP TAXA
%
% Genus: top 25 taxa by absolute aligned consensus loading.
% Phylum: all taxa.
%% ============================================================

if lower(string(analysisLevel)) == "phylum"
    idxTop = (1:nTaxa)';
else
    rankScore = abs(load_consensus(:));
    [~, idxTop] = maxk(rankScore, min(topK, nTaxa));
    [~, ordTop] = sort(abs(load_consensus(idxTop)), 'descend');
    idxTop = idxTop(ordTop);
end

taxa_top = taxaNames(idxTop);

L_top = [ ...
    load_disagree(idxTop), ...
    load_consensus(idxTop)];

Effect_top = Effect_all(idxTop, :);
Sig_top    = Sig_all(idxTop, :);

loadingNames = ["Disagreement","Consensus"];

%% ============================================================
% BUILD COMBINED MODEL MATRIX FOR PANEL C
%% ============================================================

if isfield(ascao, 'interactions') && ~isempty(ascao.interactions)

    Mcombined = ascao.factors{1}.matrix + ...
                ascao.factors{2}.matrix + ...
                ascao.interactions{1}.matrix;

else

    warning('No ascao.interactions found. Using Factor 1 + Factor 2 only for combined score plot.');

    Mcombined = ascao.factors{1}.matrix + ...
                ascao.factors{2}.matrix;

end

%% ============================================================
% INDIVIDUAL HEATMAP EXPORTS
%% ============================================================

plotTwoColumnLoadingHeatmap( ...
    L_top, ...
    taxa_top, ...
    loadingNames, ...
    fullfile(outDir, "Schubert_ASCA_loadings_heatmap"), ...
    exportResolution);

plotMethodEffectHeatmapCompact( ...
    Effect_top, ...
    taxa_top, ...
    methodNamesEffect, ...
    Sig_top, ...
    fullfile(outDir, "Schubert_methodwise_effect_heatmap"), ...
    exportResolution);

%% ============================================================
% FINAL A-E PANEL DIRECTLY FROM DATA
%% ============================================================

createFinal_CASCA_MethodEffect_Panel( ...
    ascao.factors{1}.matrix, ...
    ascao.factors{2}.matrix, ...
    Mcombined, ...
    Mtr, ...
    G, ...
    methodOrder, ...
    groupOrder, ...
    L_top, ...
    taxa_top, ...
    Effect_top, ...
    Sig_top, ...
    methodNamesEffect, ...
    fullfile(outDir, "Schubert_CASCA_method_effect_final_panel"), ...
    analysisLevel);

fprintf('\nSchubert C-ASCA + method-effect panel saved in:\n%s\n\n', outDir);

%% ============================================================
% LOCAL FUNCTION:
% Score 1 / Score 2 individual exports
%% ============================================================

function plotScores1D_BoxSquare_NoLegend(X, classLabels, pcToPlot, classOrder, outBase, exportResolution)

    X = double(X);
    classLabels = string(classLabels(:));
    classOrder  = string(classOrder(:));

    if size(X,1) ~= numel(classLabels)
        error('Number of rows in X must match number of class labels.');
    end

    X0 = X - mean(X, 1, 'omitnan');
    X0(isnan(X0)) = 0;

    [U, S, ~] = svd(X0, 'econ');

    score = U * S;
    eigvals = diag(S).^2;
    explained = 100 * eigvals ./ sum(eigvals);

    if pcToPlot > size(score,2)
        error('Requested PC%d, but only %d PCs are available.', pcToPlot, size(score,2));
    end

    y = score(:, pcToPlot);
    pcPerc = explained(pcToPlot);

    nClass = numel(classOrder);
    cmap = zeros(nClass, 3);

    for k = 1:nClass
        cmap(k,:) = getClassColor(classOrder(k));
    end

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1 1 5.2 5.2], ...
        'InvertHardcopy', 'off', ...
        'PaperPositionMode', 'auto', ...
        'Renderer', 'painters');

    ax = axes(fig, ...
        'Units', 'normalized', ...
        'Position', [0.18 0.18 0.58 0.58]);

    hold(ax, 'on');

    dotSize = 52;
    lineHalfWidth = 0.40;
    lineWidth = 2.5;

    for k = 1:nClass

        idx = classLabels == classOrder(k);

        if ~any(idx)
            warning('Class "%s" not found in labels.', classOrder(k));
            continue;
        end

        mu = mean(y(idx), 'omitnan');

        scatter(ax, k, mu, ...
            dotSize, ...
            'o', ...
            'MarkerFaceColor', cmap(k,:), ...
            'MarkerEdgeColor', cmap(k,:), ...
            'LineWidth', 0.8);

        plot(ax, ...
            [k-lineHalfWidth, k+lineHalfWidth], ...
            [mu, mu], ...
            '-', ...
            'Color', cmap(k,:), ...
            'LineWidth', lineWidth);
    end

    yline(ax, 0, '-', ...
        'Color', [0.35 0.35 0.35], ...
        'LineWidth', 1.1);

    xlim(ax, [0.5 nClass + 0.5]);

    ax.XTick = 1:nClass;
    ax.XTickLabel = classOrder;
    xtickangle(ax, 35);

    ax.FontName = 'Arial';
    ax.FontSize = 11;
    ax.LineWidth = 1.0;
    ax.Box = 'on';
    ax.TickDir = 'out';

    grid(ax, 'on');
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.28;
    ax.GridColor = [0.75 0.75 0.75];

    ylabel(ax, sprintf('PC%d score (%.1f%% explained)', pcToPlot, pcPerc), ...
        'FontName', 'Arial', ...
        'FontSize', 12, ...
        'FontWeight', 'bold');

    xlabel(ax, '');

    yl = ylim(ax);
    dy = 0.12 * max(eps, diff(yl));
    ylim(ax, [yl(1)-dy, yl(2)+dy]);

    pbaspect(ax, [1 1 1]);

    safeExportFigure(fig, outBase, exportResolution, true);

    close(fig);
end

%% ============================================================
% LOCAL FUNCTION:
% Combined model score plot individual export
%% ============================================================

function plotCombinedModelScores_BoxSquare(M, methodLabels, groupLabels, methodOrder, groupOrder, outBase, exportResolution)

    M = double(M);
    methodLabels = string(methodLabels(:));
    groupLabels  = string(groupLabels(:));
    methodOrder  = string(methodOrder(:));
    groupOrder   = string(groupOrder(:));

    if size(M,1) ~= numel(methodLabels) || size(M,1) ~= numel(groupLabels)
        error('M, methodLabels, and groupLabels must have the same number of rows.');
    end

    Xcomb = M - mean(M, 1, 'omitnan');
    Xcomb(isnan(Xcomb)) = 0;

    [~, scoreComb, ~, ~, explainedComb] = pca(Xcomb, 'Rows', 'complete');

    cmapMethods = zeros(numel(methodOrder), 3);

    for i = 1:numel(methodOrder)
        cmapMethods(i,:) = getClassColor(methodOrder(i));
    end

    figComb = figure( ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1 1 5.2 5.2], ...
        'InvertHardcopy', 'off', ...
        'PaperPositionMode', 'auto', ...
        'Renderer', 'painters');

    axComb = axes(figComb, ...
        'Units', 'normalized', ...
        'Position', [0.18 0.18 0.58 0.58]);

    hold(axComb, 'on');

    markerSize = 90;

    for i = 1:numel(methodOrder)

        idxM = methodLabels == methodOrder(i);

        for j = 1:numel(groupOrder)

            idxG = groupLabels == groupOrder(j);
            idx  = idxM & idxG;

            if ~any(idx)
                continue;
            end

            if groupOrder(j) == "Control"
                mk = 'o';
            elseif groupOrder(j) == "Case"
                mk = '^';
            else
                mk = 's';
            end

            scatter(axComb, scoreComb(idx,1), scoreComb(idx,2), ...
                markerSize, ...
                'Marker', mk, ...
                'MarkerEdgeColor', cmapMethods(i,:), ...
                'MarkerFaceColor', cmapMethods(i,:), ...
                'MarkerFaceAlpha', 0.90, ...
                'MarkerEdgeAlpha', 1.00, ...
                'LineWidth', 1.0);
        end
    end

    xline(axComb, 0, 'k-', 'LineWidth', 0.8);
    yline(axComb, 0, 'k-', 'LineWidth', 0.8);

    xlabel(axComb, sprintf('PC1 (%.1f%%)', explainedComb(1)), ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial');

    ylabel(axComb, sprintf('PC2 (%.1f%%)', explainedComb(2)), ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial');

    box(axComb, 'on');
    grid(axComb, 'on');

    axComb.GridAlpha = 0.25;
    axComb.FontName = 'Arial';
    axComb.FontSize = 11;
    axComb.LineWidth = 1.0;
    axComb.TickDir = 'out';

    axis(axComb, 'tight');

    xl = xlim(axComb);
    yl = ylim(axComb);

    dx = 0.08 * max(eps, diff(xl));
    dy = 0.08 * max(eps, diff(yl));

    xlim(axComb, [xl(1)-dx, xl(2)+dx]);
    ylim(axComb, [yl(1)-dy, yl(2)+dy]);

    pbaspect(axComb, [1 1 1]);

    hMethod = gobjects(numel(methodOrder),1);

    for i = 1:numel(methodOrder)
        hMethod(i) = scatter(axComb, nan, nan, 80, ...
            'o', ...
            'MarkerEdgeColor', cmapMethods(i,:), ...
            'MarkerFaceColor', cmapMethods(i,:), ...
            'LineWidth', 0.9);
    end

    hGroup = gobjects(2,1);

    hGroup(1) = scatter(axComb, nan, nan, 80, ...
        'o', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', ...
        'LineWidth', 0.9);

    hGroup(2) = scatter(axComb, nan, nan, 80, ...
        '^', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', ...
        'LineWidth', 0.9);

    hSep = plot(axComb, nan, nan, ...
        'LineStyle', 'none', ...
        'Marker', 'none', ...
        'Color', 'none');

    hLegend = [hMethod; hSep; hGroup(:)];

    legendLabels = [ ...
        cellstr(methodOrder(:)); ...
        {''}; ...
        {'Control'; 'Case'}];

    leg = legend(axComb, hLegend, legendLabels, ...
        'Units', 'normalized', ...
        'Location', 'eastoutside', ...
        'FontSize', 10.5, ...
        'FontName', 'Arial', ...
        'Box', 'on');

    leg.EdgeColor = [0 0 0];
    leg.LineWidth = 0.6;
    leg.ItemTokenSize = [14 12];
    title(leg, '');

    safeExportFigure(figComb, outBase, exportResolution, true);

    close(figComb);
end

%% ============================================================
% INDIVIDUAL LOADING HEATMAP
%% ============================================================

function plotTwoColumnLoadingHeatmap(L_plot, taxaNames, colNames, outBase, exportResolution)

    taxaNames = string(taxaNames(:));
    colNames  = string(colNames(:))';

    [nRows, nCols] = size(L_plot);

    showValues = true;
    nDec = 2;

    cellW = 0.60;
    cellH = 0.24;   % slightly taller than before

    if nRows <= 6
        cellH = 0.45;
    end

    leftMargin   = 2.40;
    rightMargin  = 0.55;
    bottomMargin = 0.95;
    topMargin    = 0.20;

    cbGap = 0.12;
    cbW   = 0.18;

    axW  = nCols * cellW;
    axH  = nRows * cellH;
    figW = leftMargin + axW + cbGap + cbW + rightMargin;
    figH = bottomMargin + axH + topMargin;

    figLoad = figure('Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1 1 figW figH], ...
        'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 figW figH], ...
        'PaperSize', [figW figH], ...
        'PaperPositionMode', 'auto', ...
        'InvertHardcopy', 'off', ...
        'Renderer', 'painters');

    axLoad = axes(figLoad, ...
        'Units', 'inches', ...
        'Position', [leftMargin bottomMargin axW axH]);

    imagesc(axLoad, L_plot);
    clim(axLoad, [-1 1]);
    colormap(axLoad, blueWhiteRedCMap(256));

    axLoad.XTick = 1:nCols;
    axLoad.XTickLabel = cellstr(colNames);
    axLoad.YTick = [];
    axLoad.FontName = 'Arial';
    axLoad.FontSize = 10;
    axLoad.TickLength = [0 0];
    axLoad.Box = 'off';
    axLoad.LineWidth = 0.8;

    xtickangle(axLoad, 35);

    xlim(axLoad, [0.5 nCols + 0.5]);
    ylim(axLoad, [0.5 nRows + 0.5]);
    set(axLoad, 'YDir', 'normal');

    hold(axLoad, 'on');

    drawHeatmapGrid(axLoad, nRows, nCols);

    xLabelPos = -3.25;

    for i = 1:nRows

        taxonLabel = char(taxaNames(i));

        if startsWith(string(taxonLabel), "unknown", 'IgnoreCase', true)
            fontAngle = 'normal';
        else
            fontAngle = 'italic';
        end

        text(axLoad, xLabelPos, i, taxonLabel, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontAngle', fontAngle, ...
            'FontSize', 8.8, ...
            'Color', [0 0 0], ...
            'Clipping', 'off');
    end

    if showValues
        addHeatmapNumbers(axLoad, L_plot, nDec, 0.55, 7.6);
    end

    cbLoad = colorbar(axLoad);
    cbLoad.Units = 'inches';
    cbLoad.Position = [leftMargin + axW + cbGap, bottomMargin, cbW, axH];

    cbLoad.FontName = 'Arial';
    cbLoad.FontSize = 8.8;
    cbLoad.LineWidth = 0.6;
    cbLoad.Box = 'off';
    cbLoad.Label.String = 'CASCA loading';
    cbLoad.Label.FontName = 'Arial';
    cbLoad.Label.FontSize = 10;
    cbLoad.Label.FontWeight = 'normal';
    cbLoad.Ticks = [-1 0 1];
    cbLoad.TickLabels = {'-1','0','1'};
    cbLoad.TickDirection = 'out';

    safeExportFigure(figLoad, outBase, exportResolution, true);

    close(figLoad);
end

%% ============================================================
% INDIVIDUAL METHOD-WISE EFFECT HEATMAP
%% ============================================================

function plotMethodEffectHeatmapCompact(Effect_top, taxaNames, methodNames, Sig_top, outBase, exportResolution)

    taxaNames = string(taxaNames(:));
    methodNames = string(methodNames(:))';

    [nRows, nCols] = size(Effect_top);

    nDec = 2;
    cellW = 0.60;
    cellH = 0.24;   % same row height as loading heatmap

    if nRows <= 6
        cellH = 0.45;
    end

    leftMargin   = 2.40;
    rightMargin  = 0.55;
    bottomMargin = 0.95;
    topMargin    = 0.20;

    cbGap = 0.12;
    cbW   = 0.18;

    axW  = nCols * cellW;
    axH  = nRows * cellH;
    figW = leftMargin + axW + cbGap + cbW + rightMargin;
    figH = bottomMargin + axH + topMargin;

    figEff = figure('Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1 1 figW figH], ...
        'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 figW figH], ...
        'PaperSize', [figW figH], ...
        'PaperPositionMode', 'auto', ...
        'InvertHardcopy', 'off', ...
        'Renderer', 'painters');

    axEff = axes(figEff, ...
        'Units', 'inches', ...
        'Position', [leftMargin bottomMargin axW axH]);

    imagesc(axEff, Effect_top);
    clim(axEff, [-1 1]);
    colormap(axEff, blueWhiteRedCMap(256));

    axEff.XTick = 1:nCols;
    axEff.XTickLabel = cellstr(methodNames);
    axEff.YTick = [];
    axEff.FontName = 'Arial';
    axEff.FontSize = 10;
    axEff.TickLength = [0 0];
    axEff.Box = 'off';
    axEff.LineWidth = 0.8;

    xtickangle(axEff, 35);

    xlim(axEff, [0.5 nCols + 0.5]);
    ylim(axEff, [0.5 nRows + 0.5]);
    set(axEff, 'YDir', 'normal');

    hold(axEff, 'on');

    drawHeatmapGrid(axEff, nRows, nCols);

    xLabelPos = -2.15;

    for i = 1:nRows

        taxonLabel = char(taxaNames(i));

        if startsWith(string(taxonLabel), "unknown", 'IgnoreCase', true)
            fontAngle = 'normal';
        else
            fontAngle = 'italic';
        end

        text(axEff, xLabelPos, i, taxonLabel, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontAngle', fontAngle, ...
            'FontSize', 8.8, ...
            'Color', [0 0 0], ...
            'Clipping', 'off');
    end

    addHeatmapNumbers(axEff, Effect_top, nDec, 0.55, 6.7);

    for i = 1:nRows
        for j = 1:nCols
            if Sig_top(i,j)
                text(axEff, j, i - 0.34, '*', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontName', 'Arial', ...
                    'FontWeight', 'bold', ...
                    'FontSize', 8.0, ...
                    'Color', [0 0 0], ...
                    'Clipping', 'on');
            end
        end
    end

    cbEff = colorbar(axEff);
    cbEff.Units = 'inches';
    cbEff.Position = [leftMargin + axW + cbGap, bottomMargin, cbW, axH];

    cbEff.FontName = 'Arial';
    cbEff.FontSize = 8.8;
    cbEff.LineWidth = 0.6;
    cbEff.Box = 'off';
    cbEff.Label.String = 'Median(Case) - Median(Control)';
    cbEff.Label.FontName = 'Arial';
    cbEff.Label.FontSize = 10;
    cbEff.Label.FontWeight = 'normal';
    cbEff.Ticks = [-1 0 1];
    cbEff.TickLabels = {'-1','0','1'};
    cbEff.TickDirection = 'out';

    safeExportFigure(figEff, outBase, exportResolution, true);

    close(figEff);
end

%% ============================================================
% FINAL PUBLICATION PANEL FUNCTION
%% ============================================================
function createFinal_CASCA_MethodEffect_Panel( ...
    X_factor1, X_factor2, Mcombined, ...
    methodLabels, groupLabels, methodOrder, groupOrder, ...
    L_top, taxa_top, Effect_top, Sig_top, effectMethodNames, outBase, analysisLevel)

    panelWidthIn  = 7.2;
    panelHeightIn = 8.6;
    pngDPI = 600;

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [0.25 0.25 panelWidthIn panelHeightIn], ...
        'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 panelWidthIn panelHeightIn], ...
        'PaperSize', [panelWidthIn panelHeightIn], ...
        'InvertHardcopy', 'off', ...
        'Renderer', 'painters');

    if lower(string(analysisLevel)) == "phylum"

        topY = 0.720;
        topH = 0.110;

        botY = 0.245;
        botH = 0.365;

    else

        % Move bottom row upward slightly so x-labels are not cropped.
        % Keep D and E identical in height.
        topY = 0.875;
        topH = 0.100;

        botY = 0.070;
        botH = 0.720;

    end

    posA = [0.055  topY  0.245  topH];
    posB = [0.350  topY  0.210  topH];
    posC = [0.620  topY  0.325  topH];

    posD_taxa = [0.025  botY  0.135  botH];
    posD_main = [0.170  botY  0.180  botH];
    posD_cb   = [0.364  botY  0.018  botH];

    posE_main = [0.455  botY  0.420  botH];
    posE_cb   = [0.892  botY  0.018  botH];

    axA = axes(fig, 'Units', 'normalized', 'Position', posA);
    plotPanelScore1D(axA, X_factor1, methodLabels, 1, methodOrder);
    addPanelLabel(fig, posA, 'A');

    axB = axes(fig, 'Units', 'normalized', 'Position', posB);
    plotPanelScore1D(axB, X_factor2, groupLabels, 1, groupOrder);
    addPanelLabel(fig, posB, 'B');

    axC = axes(fig, 'Units', 'normalized', 'Position', posC);
    plotPanelCombinedScores(axC, Mcombined, methodLabels, groupLabels, methodOrder, groupOrder);
    addPanelLabel(fig, posC, 'C');

    axDtaxa = axes(fig, 'Units', 'normalized', 'Position', posD_taxa);
    axDmain = axes(fig, 'Units', 'normalized', 'Position', posD_main);
    axDcb   = axes(fig, 'Units', 'normalized', 'Position', posD_cb);

    plotPanelLoadingHeatmap( ...
        axDmain, ...
        axDtaxa, ...
        axDcb, ...
        L_top, ...
        taxa_top, ...
        ["Disagreement","Consensus"]);

    addPanelLabel(fig, [posD_taxa(1), posD_main(2), posD_taxa(3)+posD_main(3), posD_main(4)], 'D');

    axEmain = axes(fig, 'Units', 'normalized', 'Position', posE_main);
    axEcb   = axes(fig, 'Units', 'normalized', 'Position', posE_cb);

    plotPanelMethodEffectHeatmap_NoTaxa( ...
        axEmain, ...
        axEcb, ...
        Effect_top, ...
        Sig_top, ...
        effectMethodNames);

    addPanelLabel(fig, posE_main, 'E');

    safeExportFigure(fig, outBase, pngDPI, true);

    close(fig);
end










%% ============================================================
% PANEL LABEL
%% ============================================================

function addPanelLabel(fig, pos, txt)

    labelW = 0.040;
    labelH = 0.022;

    x = pos(1) - 0.012;
    y = pos(2) + pos(4) + 0.004;

    if x < 0.002
        x = 0.002;
    end

    if y + labelH > 0.995
        y = 0.995 - labelH;
    end

    annotation(fig, 'textbox', ...
        [x, y, labelW, labelH], ...
        'String', txt, ...
        'LineStyle', 'none', ...
        'FontName', 'Arial', ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'Color', 'k', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'FitBoxToText', 'off');
end

%% ============================================================
% PANEL A / B: 1D SCORE PLOT
%% ============================================================

function plotPanelScore1D(ax, X, classLabels, pcToPlot, classOrder)

    X = double(X);
    classLabels = string(classLabels(:));
    classOrder  = string(classOrder(:));

    X0 = X - mean(X, 1, 'omitnan');
    X0(isnan(X0)) = 0;

    [U, S, ~] = svd(X0, 'econ');

    score = U * S;
    eigvals = diag(S).^2;
    explained = 100 * eigvals ./ sum(eigvals);

    y = score(:, pcToPlot);
    pcPerc = explained(pcToPlot);

    hold(ax, 'on');

    nClass = numel(classOrder);
    cmap = zeros(nClass, 3);

    for k = 1:nClass
        cmap(k,:) = getClassColor(classOrder(k));
    end

    dotSize = 30;
    lineHalfWidth = 0.28;
    lineWidth = 1.8;

    for k = 1:nClass

        idx = classLabels == classOrder(k);

        if ~any(idx)
            continue;
        end

        mu = mean(y(idx), 'omitnan');

        scatter(ax, k, mu, ...
            dotSize, ...
            'o', ...
            'MarkerFaceColor', cmap(k,:), ...
            'MarkerEdgeColor', cmap(k,:), ...
            'LineWidth', 0.7);

        plot(ax, ...
            [k-lineHalfWidth, k+lineHalfWidth], ...
            [mu, mu], ...
            '-', ...
            'Color', cmap(k,:), ...
            'LineWidth', lineWidth);
    end

    yline(ax, 0, '-', ...
        'Color', [0.35 0.35 0.35], ...
        'LineWidth', 0.8);

    xlim(ax, [0.5 nClass + 0.5]);

    ax.XTick = 1:nClass;
    ax.XTickLabel = cellstr(classOrder);
    xtickangle(ax, 35);

    ax.FontName = 'Arial';
    ax.FontSize = 7.4;
    ax.LineWidth = 0.75;
    ax.Box = 'on';
    ax.TickDir = 'out';
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.25;
    ax.GridColor = [0.75 0.75 0.75];

    ylabel(ax, sprintf('PC%d score\n(%.1f%% explained)', pcToPlot, pcPerc), ...
        'FontName', 'Arial', ...
        'FontSize', 8.0, ...
        'FontWeight', 'bold');

    xlabel(ax, '');

    yl = ylim(ax);
    dy = 0.10 * max(eps, diff(yl));
    ylim(ax, [yl(1)-dy, yl(2)+dy]);

    pbaspect(ax, [1 0.92 1]);
end

%% ============================================================
% PANEL C: COMBINED MODEL SCORE PLOT
%% ============================================================

function plotPanelCombinedScores(ax, M, methodLabels, groupLabels, methodOrder, groupOrder)

    M = double(M);
    methodLabels = string(methodLabels(:));
    groupLabels  = string(groupLabels(:));
    methodOrder  = string(methodOrder(:));
    groupOrder   = string(groupOrder(:));

    Xcomb = M - mean(M, 1, 'omitnan');
    Xcomb(isnan(Xcomb)) = 0;

    [~, scoreComb, ~, ~, explainedComb] = pca(Xcomb, 'Rows', 'complete');

    hold(ax, 'on');

    cmapMethods = zeros(numel(methodOrder), 3);

    for i = 1:numel(methodOrder)
        cmapMethods(i,:) = getClassColor(methodOrder(i));
    end

    markerSize = 42;

    for i = 1:numel(methodOrder)

        idxM = methodLabels == methodOrder(i);

        for j = 1:numel(groupOrder)

            idxG = groupLabels == groupOrder(j);
            idx  = idxM & idxG;

            if ~any(idx)
                continue;
            end

            if groupOrder(j) == "Control"
                mk = 'o';
            else
                mk = '^';
            end

            scatter(ax, scoreComb(idx,1), scoreComb(idx,2), ...
                markerSize, ...
                'Marker', mk, ...
                'MarkerEdgeColor', cmapMethods(i,:), ...
                'MarkerFaceColor', cmapMethods(i,:), ...
                'MarkerFaceAlpha', 0.92, ...
                'MarkerEdgeAlpha', 1.00, ...
                'LineWidth', 0.7);
        end
    end

    xline(ax, 0, 'k-', 'LineWidth', 0.65);
    yline(ax, 0, 'k-', 'LineWidth', 0.65);

    xlabel(ax, sprintf('PC1 (%.1f%%)', explainedComb(1)), ...
        'FontSize', 8.5, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial');

    ylabel(ax, sprintf('PC2 (%.1f%%)', explainedComb(2)), ...
        'FontSize', 8.5, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial');

    box(ax, 'on');
    grid(ax, 'on');

    ax.GridAlpha = 0.25;
    ax.FontName = 'Arial';
    ax.FontSize = 7.5;
    ax.LineWidth = 0.75;
    ax.TickDir = 'out';

    axis(ax, 'tight');

    xl = xlim(ax);
    yl = ylim(ax);

    dx = 0.08 * max(eps, diff(xl));
    dy = 0.08 * max(eps, diff(yl));

    xlim(ax, [xl(1)-dx, xl(2)+dx]);
    ylim(ax, [yl(1)-dy, yl(2)+dy]);

    pbaspect(ax, [1 1 1]);

    hMethod = gobjects(numel(methodOrder),1);

    for i = 1:numel(methodOrder)
        hMethod(i) = scatter(ax, nan, nan, 32, ...
            'o', ...
            'MarkerEdgeColor', cmapMethods(i,:), ...
            'MarkerFaceColor', cmapMethods(i,:), ...
            'LineWidth', 0.7);
    end

    hGroup = gobjects(2,1);

    hGroup(1) = scatter(ax, nan, nan, 32, ...
        'o', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', ...
        'LineWidth', 0.7);

    hGroup(2) = scatter(ax, nan, nan, 32, ...
        '^', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', ...
        'LineWidth', 0.7);

    hSep = plot(ax, nan, nan, ...
        'LineStyle', 'none', ...
        'Marker', 'none', ...
        'Color', 'none');

    hLegend = [hMethod; hSep; hGroup(:)];

    legendLabels = [ ...
        cellstr(methodOrder(:)); ...
        {''}; ...
        {'Control'; 'Case'}];

    leg = legend(ax, hLegend, legendLabels, ...
        'Location', 'eastoutside', ...
        'FontSize', 7, ...
        'FontName', 'Arial', ...
        'Box', 'on');

    leg.EdgeColor = [0 0 0];
    leg.LineWidth = 0.4;
    leg.ItemTokenSize = [10 8];
end

%% ============================================================
% PANEL D: C-ASCA LOADING HEATMAP
%% ============================================================

function plotPanelLoadingHeatmap(axMain, axTaxa, axCB, L_plot, taxaNames, colNames)

    taxaNames = string(taxaNames(:));
    colNames = string(colNames(:))';

    nRows = size(L_plot,1);
    nCols = size(L_plot,2);

    imagesc(axMain, L_plot);
    clim(axMain, [-1 1]);
    colormap(axMain, blueWhiteRedCMap(256));

    axMain.XTick = 1:nCols;
    axMain.XTickLabel = cellstr(colNames);
    axMain.YTick = [];
    xtickangle(axMain, 35);

    axMain.FontName = 'Arial';
    axMain.FontSize = 8;
    axMain.TickLength = [0 0];
    axMain.Box = 'off';
    axMain.LineWidth = 0.6;

    xlim(axMain, [0.5 nCols + 0.5]);
    ylim(axMain, [0.5 nRows + 0.5]);
    set(axMain, 'YDir', 'normal');

    hold(axMain, 'on');

    drawHeatmapGrid(axMain, nRows, nCols);
    addHeatmapNumbers(axMain, L_plot, 2, 0.55, 7.2);

    cla(axTaxa);
    axis(axTaxa, 'off');
    xlim(axTaxa, [0 1]);
    ylim(axTaxa, [0.5 nRows + 0.5]);
    set(axTaxa, 'YDir', 'normal');

    for i = 1:nRows

        taxonLabel = char(taxaNames(i));

        if startsWith(string(taxonLabel), "unknown", 'IgnoreCase', true)
            fontAngle = 'normal';
        else
            fontAngle = 'italic';
        end

        text(axTaxa, 0.00, i, taxonLabel, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontAngle', fontAngle, ...
            'FontSize', 7.6, ...
            'Color', [0 0 0], ...
            'Clipping', 'off');
    end

    cb = colorbar(axMain);
    cb.Units = 'normalized';
    cb.Position = axCB.Position;

    cb.FontName = 'Arial';
    cb.FontSize = 8;
    cb.LineWidth = 0.5;
    cb.Box = 'off';
    cb.Ticks = [-1 0 1];
    cb.TickLabels = {'-1','0','1'};
    cb.Label.String = 'CASCA loading';
    cb.Label.FontName = 'Arial';
    cb.Label.FontSize = 9;
    cb.Label.FontWeight = 'normal';

    axis(axCB, 'off');
end

%% ============================================================
% PANEL E: METHOD-WISE EFFECT HEATMAP WITHOUT TAXA LABELS
%% ============================================================
function plotPanelMethodEffectHeatmap_NoTaxa(axMain, axCB, Effect_top, Sig_top, methodNames)

    methodNames = string(methodNames(:))';

    nRows = size(Effect_top,1);
    nCols = size(Effect_top,2);

    imagesc(axMain, Effect_top);
    clim(axMain, [-1 1]);
    colormap(axMain, blueWhiteRedCMap(256));

    axMain.XTick = 1:nCols;
    axMain.XTickLabel = cellstr(methodNames);
    axMain.YTick = [];
    xtickangle(axMain, 35);

    axMain.FontName = 'Arial';
    axMain.FontSize = 7.8;
    axMain.TickLength = [0 0];
    axMain.Box = 'off';
    axMain.LineWidth = 0.6;

    xlim(axMain, [0.5 nCols + 0.5]);
    ylim(axMain, [0.5 nRows + 0.5]);
    set(axMain, 'YDir', 'normal');

    hold(axMain, 'on');

    drawHeatmapGrid(axMain, nRows, nCols);

    % Numeric values centered in each cell
    addHeatmapNumbers(axMain, Effect_top, 2, 0.55, 6.3);

    % Put significance star ABOVE the numeric value, not below it
    for i = 1:nRows
        for j = 1:nCols
            if Sig_top(i,j)
                text(axMain, j, i + 0.22, '*', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontName', 'Arial', ...
                    'FontWeight', 'bold', ...
                    'FontSize', 7.2, ...
                    'Color', [0 0 0], ...
                    'Clipping', 'on');
            end
        end
    end

    cb = colorbar(axMain);
    cb.Units = 'normalized';
    cb.Position = axCB.Position;

    cb.FontName = 'Arial';
    cb.FontSize = 8;
    cb.LineWidth = 0.5;
    cb.Box = 'off';
    cb.Ticks = [-1 0 1];
    cb.TickLabels = {'-1','0','1'};
    cb.Label.String = 'Median(Case) - Median(Control)';
    cb.Label.FontName = 'Arial';
    cb.Label.FontSize = 9;
    cb.Label.FontWeight = 'normal';

    axis(axCB, 'off');
end
%% ============================================================
% LOCAL FUNCTION:
% Kruskal-Wallis + BH-FDR
%% ============================================================

function [pvals, qvals, isSig, effect] = local_kw_bh_feature_test(Xblock, g, alpha)

    if nargin < 3
        alpha = 0.05;
    end

    g = g(:);

    if size(Xblock,1) ~= numel(g)
        error('Rows of Xblock must match group vector length.');
    end

    p = size(Xblock,2);

    pvals = nan(p,1);
    qvals = nan(p,1);
    isSig = false(p,1);
    effect = nan(p,1);

    for j = 1:p

        x = Xblock(:,j);

        valid = ~isnan(x) & ~isnan(g);
        xj = x(valid);
        gj = g(valid);

        if numel(unique(gj)) < 2
            pvals(j) = NaN;
            effect(j) = NaN;
            continue;
        end

        if numel(unique(xj)) <= 1
            pvals(j) = 1;
            effect(j) = 0;
            continue;
        end

        pvals(j) = kruskalwallis(xj, gj, 'off');

        x_case = xj(gj == 1);
        x_ctrl = xj(gj == -1);

        effect(j) = median(x_case, 'omitnan') - median(x_ctrl, 'omitnan');
    end

    qvals = local_bh_fdr(pvals);
    isSig = qvals <= alpha;
end

%% ============================================================
% LOCAL FUNCTION:
% BH-FDR
%% ============================================================

function q = local_bh_fdr(p)

    p = p(:);
    q = nan(size(p));

    valid = ~isnan(p);
    pv = p(valid);

    if isempty(pv)
        return;
    end

    [pv_sorted, order] = sort(pv, 'ascend');
    m = numel(pv_sorted);

    q_sorted = pv_sorted .* m ./ (1:m)';

    for i = m-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end

    q_sorted(q_sorted > 1) = 1;

    q_valid = nan(size(pv));
    q_valid(order) = q_sorted;

    q(valid) = q_valid;
end

%% ============================================================
% LOCAL FUNCTION:
% Draw heatmap grid
%% ============================================================

function drawHeatmapGrid(ax, nRows, nCols)

    for r = 0.5:1:(nRows + 0.5)
        plot(ax, [0.5 nCols + 0.5], [r r], '-', ...
            'Color', [0.88 0.88 0.88], ...
            'LineWidth', 0.25);
    end

    for c = 0.5:1:(nCols + 0.5)
        plot(ax, [c c], [0.5 nRows + 0.5], '-', ...
            'Color', [0.88 0.88 0.88], ...
            'LineWidth', 0.25);
    end
end

%% ============================================================
% LOCAL FUNCTION:
% Add heatmap numbers
%% ============================================================

function addHeatmapNumbers(ax, M, nDec, whiteThreshold, fontSize)

    [nRows, nCols] = size(M);

    for i = 1:nRows
        for j = 1:nCols

            val = M(i,j);

            if abs(val) < 5e-4
                labelStr = '0';
            else
                labelStr = sprintf(['%.' num2str(nDec) 'f'], val);
            end

            if abs(val) > whiteThreshold
                txtColor = [1 1 1];
            else
                txtColor = [0 0 0];
            end

            text(ax, j, i, labelStr, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontName', 'Arial', ...
                'FontSize', fontSize, ...
                'Color', txtColor);
        end
    end
end

%% ============================================================
% LOCAL FUNCTION:
% Safe export
%% ============================================================

function safeExportFigure(figHandle, outBase, dpi, makePdf)

    if nargin < 4
        makePdf = true;
    end

    outBase = string(outBase);
    outFolder = fileparts(char(outBase));

    if ~isempty(outFolder) && ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    pngFile = char(outBase + ".png");
    pdfFile = char(outBase + ".pdf");

    drawnow;

    exportgraphics(figHandle, pngFile, ...
        'Resolution', dpi, ...
        'BackgroundColor', 'white');

    if makePdf
        exportgraphics(figHandle, pdfFile, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    end
end

%% ============================================================
% LOCAL FUNCTION:
% Class colors
%% ============================================================

function rgb = getClassColor(className)

    className = string(className);

    switch className

        case "Raw Data"
            rgb = [0.0000 0.6000 0.2500];

        case "TSS"
            rgb = [0.6350 0.0780 0.1840];

        case "Rarefaction"
            rgb = [0.9290 0.6940 0.1250];

        case "CLR"
            rgb = [0.4940 0.1840 0.5560];

        case "Control"
            rgb = [0.0000 0.4470 0.7410];

        case "Case"
            rgb = [0.8500 0.3250 0.0980];

        otherwise
            rgb = [0.3000 0.3000 0.3000];

    end
end

%% ============================================================
% LOCAL FUNCTION:
% Blue-white-red colormap
%% ============================================================

function cmap = blueWhiteRedCMap(n)

    if nargin < 1
        n = 256;
    end

    half = floor(n/2);

    blue  = [0.000 0.447 0.741];
    white = [1.000 1.000 1.000];
    red   = [0.850 0.000 0.000];

    cmap1 = [linspace(blue(1),  white(1), half)', ...
             linspace(blue(2),  white(2), half)', ...
             linspace(blue(3),  white(3), half)'];

    cmap2 = [linspace(white(1), red(1), n-half)', ...
             linspace(white(2), red(2), n-half)', ...
             linspace(white(3), red(3), n-half)'];

    cmap = [cmap1; cmap2];
end