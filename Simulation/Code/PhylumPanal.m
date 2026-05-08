%% ============================================================
%  FULL C-ASCA + oMEDA SCRIPT + FINAL A-E PUBLICATION PANEL
%  PHYLUM VERSION
%
%  Final panel:
%    A = Factor 1 / method-related score plot
%    B = Factor 2 / group-related score plot
%    C = Combined model score plot
%    D = C-ASCA loading heatmap WITH taxa labels
%        colour = signed loading value
%        number = absolute loading value
%        star   = parglmMC significance
%
%    E = oMEDA over PCA heatmap WITHOUT taxa labels
%        colour = signed oMEDA value
%        number = absolute oMEDA value
%
%  Important:
%    - oMEDA normalization is kept exactly as before using normalizeVector.
%    - Heatmap colours remain signed.
%    - Cell numbers are absolute values only.
%    - 0.00 and -0.00 are printed as 0.
%    - Panel D stars come from parglmMC.
%
%  Export:
%    - individual plots: PNG 1200 dpi + vector PDF
%    - final panel: PNG 1200 dpi + vector PDF
%% ============================================================

clearvars;
close all;
clc;

%% ============================================================
% USER OPTION
%% ============================================================

analysisLevel = "phylum";   % fixed for this script

%% ============================================================
% LOAD DATA AND SET OUTPUT FOLDER
%% ============================================================

switch lower(analysisLevel)

    case "genus"
        load genus.mat;
        load GenusASCA.mat;
        outDir = fullfile(pwd, 'Final_Genus_CASCA_oMEDA_Figures');

    case "phylum"
        load phylum.mat;
        load PhylumASCA.mat;
        outDir = fullfile(pwd, 'Final_Phylum_CASCA_oMEDA_Figures');

    otherwise
        error('analysisLevel must be "genus" or "phylum".');
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

exportResolution = 1200;

%% ============================================================
% BUILD ASCA INPUT
%% ============================================================

X = table2array(Datos(:, 1:(size(Datos,2)-3)));

depth = sum(X, 2);

Method = string(Datos.Methods);

Method(Method == "raw data")    = "Raw Data";
Method(Method == "Raw Data")    = "Raw Data";
Method(Method == "TSS")         = "TSS";
Method(Method == "rarefaction") = "Rarefaction";
Method(Method == "Rarefaction") = "Rarefaction";
Method(Method == "CLR")         = "CLR";
Method(Method == "rCLR")        = "CLR";

MethodCode = zeros(size(Method));
MethodCode(Method == "Raw Data")    = 1;
MethodCode(Method == "TSS")         = 2;
MethodCode(Method == "Rarefaction") = 3;
MethodCode(Method == "CLR")         = 4;

if any(MethodCode == 0)
    error('Some methods were not recognized. Check Datos.Methods.');
end

Groups = Datos.Groups;

Group = ones(size(Groups));

if isnumeric(Groups)

    Group(Groups == 1) = -1;
    Group(Groups == 2) =  1;

else

    Groups = string(Groups);

    Group(Groups == "Healthy")   = -1;
    Group(Groups == "Control")   = -1;
    Group(Groups == "Case")      =  1;
    Group(Groups == "sick")      =  1;
    Group(Groups == "Dysbiosis") =  1;

end

observation = Datos.observation;

%% ============================================================
% TAXA LABELS
%% ============================================================

% Correct phylum names
lab = {
    'Bacillota'
    'Bacteroidota'
    'Actinomycetota'
    'Verrucomicrobiota'
    'Pseudomonadota'
};

if exist('lab', 'var')
    taxaNames = string(lab(:));
else
    taxaNames = string(Datos.Properties.VariableNames(1:(size(Datos,2)-3)));
end

nTaxa = numel(taxaNames);

if size(X,2) ~= nTaxa
    error('Number of taxa labels does not match number of columns in X.');
end

%% ============================================================
% BLOCK-WISE SSQ NORMALIZATION
%% ============================================================

Y = X;
num_parts = 4;
row_size = size(X,1) / num_parts;
d = depth;

if mod(size(X,1), num_parts) ~= 0
    error('Number of rows in X must be divisible by 4.');
end

for i = 1:num_parts

    idx = (i-1)*row_size + 1 : i*row_size;

    part = X(idx, :);
    part2 = preprocess2D(part, 'Preprocessing', 1);

    sum_sq = sum(part2(:).^2);

    if sum_sq > 0
        part = part / sqrt(sum_sq);
    end

    Y(idx, :) = part;
    d(idx, :) = tiedrank(depth(idx, :));

end

X = Y;
depth = d;

%% ============================================================
% RUN ASCA / PARGLM
%% ============================================================

F = [MethodCode, Group, observation, depth];

[T, struct] = parglm(X, F, ...
    'Preprocessing', 1, ...
    'Model', [1 2], ...
    'Ordinal', [0 0 0 1], ...
    'Random', [0 0 1 0], ...
    'Nested', [2,3]);

disp("ASCA with SSQ completed.");

if height(T) >= 4
    T.Source(1:4) = {'A', 'B', 'C(B)', 'D'};
end

ascao = asca(struct);
model = ascao;

%% ============================================================
% parglmMC SIGNIFICANCE FOR PANEL D
%
% Factor 1 = Method / disagreement
% Factor 2 = Group / consensus
%
% These p-values are used only as stars in the C-ASCA loading heatmap.
%% ============================================================

alpha = 0.01;
nPermMC = 1000;   % use 10 only for fast testing; use 1000 for final figures

[T_mc, parglmMCStruct] = parglmMC(X, F, ...
    'Preprocessing', 1, ...
    'Model', [1 2], ...
    'Ordinal', [0 0 0 1], ...
    'Random',  [0 0 1 0], ...
    'Nested',  [2 3], ...
    'Permutations', nPermMC, ...
    'Mtc', 3);

disp('parglmMC table:');
disp(T_mc);

if ~isfield(parglmMCStruct, 'p')
    error('parglmMCStruct does not contain field p.');
end

if size(parglmMCStruct.p,1) ~= nTaxa
    error('parglmMCStruct.p rows must match number of taxa.');
end

if size(parglmMCStruct.p,2) < 2
    error('parglmMCStruct.p must contain at least two columns: Method and Group.');
end

p_disagreement = parglmMCStruct.p(:,1);   % Method / disagreement
p_consensus    = parglmMCStruct.p(:,2);   % Group / consensus

SigLoad = [ ...
    p_disagreement(:) <= alpha, ...
    p_consensus(:)    <= alpha];

fprintf('\nparglmMC p-value diagnostics\n');

fprintf('Disagreement p: min %.3g | median %.3g | max %.3g | n<=%.3f %d/%d\n', ...
    min(p_disagreement,[],'omitnan'), ...
    median(p_disagreement,'omitnan'), ...
    max(p_disagreement,[],'omitnan'), ...
    alpha, ...
    sum(p_disagreement <= alpha), ...
    numel(p_disagreement));

fprintf('Consensus p:    min %.3g | median %.3g | max %.3g | n<=%.3f %d/%d\n', ...
    min(p_consensus,[],'omitnan'), ...
    median(p_consensus,'omitnan'), ...
    max(p_consensus,[],'omitnan'), ...
    alpha, ...
    sum(p_consensus <= alpha), ...
    numel(p_consensus));

%% ============================================================
% LABELS FOR PLOTTING
%% ============================================================

G = string(Group);
G(G == "-1") = "Control";
G(G == "1")  = "Case";

Mtr = Method;

methodOrder = ["Raw Data","TSS","Rarefaction","CLR"];
groupOrder  = ["Control","Case"];

%% ============================================================
% SCORE PLOT 1: FACTOR 1 / DISAGREEMENT
%% ============================================================

plotScores1D_BoxSquare_NoLegend( ...
    ascao.factors{1}.matrix, ...
    Mtr, ...
    1, ...
    methodOrder, ...
    fullfile(outDir, "Factor1_disagreement_scores_by_method"), ...
    exportResolution);

%% ============================================================
% SCORE PLOT 2: FACTOR 2 / CONSENSUS
%% ============================================================

plotScores1D_BoxSquare_NoLegend( ...
    ascao.factors{2}.matrix, ...
    G, ...
    1, ...
    groupOrder, ...
    fullfile(outDir, "Factor2_consensus_scores_by_group"), ...
    exportResolution);

%% ============================================================
% SCORE PLOT 3: COMBINED MODEL
%% ============================================================

if ~isfield(ascao, 'interactions') || isempty(ascao.interactions)

    error('No interaction matrix found in ascao.interactions. Cannot create combined model score plot.');

else

    Mcombined = ascao.factors{1}.matrix + ...
                ascao.factors{2}.matrix + ...
                ascao.interactions{1}.matrix;

    plotCombinedModelScores_BoxSquare( ...
        Mcombined, ...
        Mtr, ...
        G, ...
        methodOrder, ...
        groupOrder, ...
        fullfile(outDir, "combined_model_scores"), ...
        exportResolution);
end

%% ============================================================
% LOADING HEATMAP: C-ASCA LOADINGS
%% ============================================================

load_disagree_raw  = ascao.factors{1}.loads(:,1);
load_consensus_raw = ascao.factors{2}.loads(:,1);

if numel(taxaNames) ~= numel(load_disagree_raw) || numel(taxaNames) ~= numel(load_consensus_raw)
    error('taxaNames, load_disagree_raw, and load_consensus_raw must have the same length.');
end

% Keep disagreement as extracted
load_disagree = load_disagree_raw;

% Keep your previous consensus sign-orientation rule
a = ascao.factors{2}.scores(1);
b = Group(1);

if sign(a) == sign(b)
    load_consensus = load_consensus_raw;
else
    load_consensus = -load_consensus_raw;
end

L = [load_disagree(:), load_consensus(:)];

meanCenterDisplayedValues = true;

if meanCenterDisplayedValues
    L_plot = L - mean(L, 1, 'omitnan');
else
    L_plot = L;
end

tol = 5e-4;
L_plot(abs(L_plot) < tol) = 0;

plotTwoColumnLoadingHeatmap( ...
    L_plot, ...
    taxaNames, ...
    ["Disagreement","Consensus"], ...
    SigLoad, ...
    fullfile(outDir, "loadings_heatmap_2col"), ...
    exportResolution);

%% ============================================================
% oMEDA HEATMAP: oMEDA OVER PCA
%
% Normalization is kept exactly as before.
%% ============================================================

Nblock = row_size;

rj1 = omedaPca(X(1:Nblock,:), ...
               [1], ...
               X(1:Nblock,:), ...
               Group(1:Nblock), ...
               'Preprocessing', 1);

rj2 = omedaPca(X(Nblock+1:2*Nblock,:), ...
               [1], ...
               X(Nblock+1:2*Nblock,:), ...
               Group(1:Nblock), ...
               'Preprocessing', 1);

rj3 = omedaPca(X(2*Nblock+1:3*Nblock,:), ...
               [1], ...
               X(2*Nblock+1:3*Nblock,:), ...
               Group(1:Nblock), ...
               'Preprocessing', 1);

rj4 = omedaPca(X(3*Nblock+1:4*Nblock,:), ...
               [1], ...
               X(3*Nblock+1:4*Nblock,:), ...
               Group(1:Nblock), ...
               'Preprocessing', 1);

rj5 = omedaPca(ascao.factors{2}.matrix, ...
               [1], ...
               ascao.factors{2}.matrix, ...
               Group, ...
               'Preprocessing', 1);

rj1 = normalizeVector(rj1);
rj2 = normalizeVector(rj2);
rj3 = normalizeVector(rj3);
rj4 = normalizeVector(rj4);
rj5 = normalizeVector(rj5);

O = [rj1(:), rj2(:), rj3(:), rj4(:), rj5(:)];
O(abs(O) < tol) = 0;

if size(O,1) ~= numel(taxaNames)
    error('oMEDA matrix row count does not match taxa names.');
end

plotOmedaHeatmapCompact( ...
    O, ...
    ["Raw Data","TSS","Rarefaction","CLR","C-ASCA"], ...
    fullfile(outDir, "oMEDA_heatmap_compact"), ...
    exportResolution);

%% ============================================================
% FINAL A-E PANEL DIRECTLY FROM DATA
%% ============================================================

createFinal_CASCA_oMEDA_Panel( ...
    ascao.factors{1}.matrix, ...
    ascao.factors{2}.matrix, ...
    Mcombined, ...
    Mtr, ...
    G, ...
    methodOrder, ...
    groupOrder, ...
    L_plot, ...
    taxaNames, ...
    SigLoad, ...
    O, ...
    ["Raw Data","TSS","Rarefaction","CLR","C-ASCA"], ...
    fullfile(outDir, "CASCA_oMEDA_phylum_final_panel"), ...
    analysisLevel);

fprintf('\nAll scores, heatmaps, and final A-E panel saved in:\n%s\n\n', outDir);

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

    exportgraphics(fig, outBase + ".png", ...
        'Resolution', exportResolution, ...
        'BackgroundColor', 'white');

    exportgraphics(fig, outBase + ".pdf", ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

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

    exportgraphics(figComb, outBase + ".png", ...
        'Resolution', exportResolution, ...
        'BackgroundColor', 'white');

    exportgraphics(figComb, outBase + ".pdf", ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    close(figComb);
end

%% ============================================================
% LOCAL FUNCTION:
% Individual two-column C-ASCA loading heatmap
%% ============================================================

function plotTwoColumnLoadingHeatmap(L_plot, taxaNames, colNames, SigMask, outBase, exportResolution)

    taxaNames = string(taxaNames(:));
    colNames  = string(colNames(:))';

    [nRows, nCols] = size(L_plot);

    if nRows ~= numel(taxaNames)
        error('Number of rows in L_plot must match taxaNames.');
    end

    showValues = true;
    nDec = 2;

    cellW = 0.52;
    cellH = 0.45;

    leftMargin   = 1.65;
    rightMargin  = 0.65;
    bottomMargin = 1.10;
    topMargin    = 0.12;

    cbGap = 0.06;
    cbW   = 0.22;

    axW  = nCols * cellW;
    axH  = nRows * cellH;
    figW = leftMargin + axW + cbGap + cbW + rightMargin;
    figH = bottomMargin + axH + topMargin;

    figLoad = figure('Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1 1 figW figH], ...
        'PaperPositionMode', 'auto', ...
        'Renderer', 'painters');

    axLoad = axes(figLoad, ...
        'Units', 'inches', ...
        'Position', [leftMargin bottomMargin axW axH]);

    imagesc(axLoad, L_plot);

    clim(axLoad, [-1 1]);
    colormap(axLoad, blueWhiteRedCMap(256));

    axLoad.XTick = 1:nCols;
    axLoad.XTickLabel = cellstr(colNames);
    axLoad.YTick = 1:nRows;
    axLoad.YTickLabel = [];

    axLoad.FontName = 'Arial';
    axLoad.FontSize = 8.5;
    axLoad.TickLength = [0 0];
    axLoad.Box = 'off';
    axLoad.LineWidth = 0.8;
    axLoad.Layer = 'top';

    xtickangle(axLoad, 35);

    xlim(axLoad, [0.5 nCols + 0.5]);
    ylim(axLoad, [0.5 nRows + 0.5]);
    set(axLoad, 'YDir', 'normal');

    hold(axLoad, 'on');

    drawHeatmapGrid(axLoad, nRows, nCols, 0.50);

    xLabelPos = -1.55;

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
        addHeatmapNumbersAbs(axLoad, L_plot, nDec, 0.55, 7.6);
    end

    addStarMarkers(axLoad, SigMask, 8.0, 0.23);

    cbLoad = colorbar(axLoad);
    cbLoad.Units = 'inches';
    cbLoad.Position = [leftMargin + axW + cbGap, bottomMargin, cbW, axH];

    cbLoad.FontName = 'Arial';
    cbLoad.FontSize = 8.5;
    cbLoad.FontWeight = 'normal';
    cbLoad.LineWidth = 0.6;
    cbLoad.Box = 'off';

    cbLoad.Label.String = 'C-ASCA loadings';
    cbLoad.Label.FontName = 'Arial';
    cbLoad.Label.FontSize = 12;
    cbLoad.Label.FontWeight = 'normal';

    cbLoad.Ticks = [-1 0 1];
    cbLoad.TickLabels = {'-1','0','1'};
    cbLoad.TickDirection = 'out';

    title(axLoad, '');
    xlabel(axLoad, '');

    exportgraphics(figLoad, outBase + ".png", ...
        'Resolution', exportResolution, ...
        'BackgroundColor', 'white');

    exportgraphics(figLoad, outBase + ".pdf", ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    close(figLoad);
end

%% ============================================================
% LOCAL FUNCTION:
% Individual compact oMEDA heatmap WITHOUT taxa labels
%% ============================================================

function plotOmedaHeatmapCompact(O, methodNames, outBase, exportResolution)

    methodNames = string(methodNames(:))';

    [nRows, nCols] = size(O);

    if nCols ~= numel(methodNames)
        error('Number of columns in O must match methodNames.');
    end

    showValues = true;
    nDec = 2;

    cellW = 0.48;
    cellH = 0.45;

    leftMargin   = 0.25;
    rightMargin  = 0.75;
    bottomMargin = 1.10;
    topMargin    = 0.12;

    cbGap = 0.06;
    cbW   = 0.22;

    axW  = nCols * cellW;
    axH  = nRows * cellH;
    figW = leftMargin + axW + cbGap + cbW + rightMargin;
    figH = bottomMargin + axH + topMargin;

    figOmeda = figure('Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1 1 figW figH], ...
        'PaperPositionMode', 'auto', ...
        'Renderer', 'painters');

    axOmeda = axes(figOmeda, ...
        'Units', 'inches', ...
        'Position', [leftMargin bottomMargin axW axH]);

    imagesc(axOmeda, O);

    clim(axOmeda, [-1 1]);
    colormap(axOmeda, blueWhiteRedCMap(256));

    axOmeda.XTick = 1:nCols;
    axOmeda.XTickLabel = cellstr(methodNames);
    axOmeda.YTick = [];
    axOmeda.YTickLabel = [];

    axOmeda.FontName = 'Arial';
    axOmeda.FontSize = 8.5;
    axOmeda.TickLength = [0 0];
    axOmeda.Box = 'off';
    axOmeda.LineWidth = 0.8;
    axOmeda.Layer = 'top';

    xtickangle(axOmeda, 35);

    xlim(axOmeda, [0.5 nCols + 0.5]);
    ylim(axOmeda, [0.5 nRows + 0.5]);
    set(axOmeda, 'YDir', 'normal');

    hold(axOmeda, 'on');

    drawHeatmapGrid(axOmeda, nRows, nCols, 0.50);

    if showValues
        addHeatmapNumbersAbs(axOmeda, O, nDec, 0.65, 7.4);
    end

    cbOmeda = colorbar(axOmeda);
    cbOmeda.Units = 'inches';
    cbOmeda.Position = [leftMargin + axW + cbGap, bottomMargin, cbW, axH];

    cbOmeda.FontName = 'Arial';
    cbOmeda.FontSize = 8.5;
    cbOmeda.FontWeight = 'normal';
    cbOmeda.LineWidth = 0.6;
    cbOmeda.Box = 'off';

    cbOmeda.Label.String = 'oMEDA over PCA';
    cbOmeda.Label.FontName = 'Arial';
    cbOmeda.Label.FontSize = 12;
    cbOmeda.Label.FontWeight = 'normal';

    cbOmeda.Ticks = [-1 0 1];
    cbOmeda.TickLabels = {'-1','0','1'};
    cbOmeda.TickDirection = 'out';

    title(axOmeda, '');
    xlabel(axOmeda, '');

    exportgraphics(figOmeda, outBase + ".png", ...
        'Resolution', exportResolution, ...
        'BackgroundColor', 'white');

    exportgraphics(figOmeda, outBase + ".pdf", ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    close(figOmeda);
end

%% ============================================================
% FINAL PUBLICATION PANEL FUNCTION
%% ============================================================

function createFinal_CASCA_oMEDA_Panel( ...
    X_factor1, X_factor2, Mcombined, ...
    methodLabels, groupLabels, methodOrder, groupOrder, ...
    L_plot, taxaNames, SigLoad, O, omedaMethodNames, outBase, analysisLevel)

    panelWidthIn  = 7.2;
    panelHeightIn = 7.8;
    pngDPI = 1200;

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

        topY = 0.725;
        topH = 0.110;

        botY = 0.255;
        botH = 0.345;

    else

        topY = 0.865;
        topH = 0.110;

        botY = 0.090;
        botH = 0.690;

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
        L_plot, ...
        taxaNames, ...
        ["Disagreement","Consensus"], ...
        SigLoad);

    addPanelLabel(fig, [posD_taxa(1), posD_main(2), posD_taxa(3)+posD_main(3), posD_main(4)], 'D');

    axEmain = axes(fig, 'Units', 'normalized', 'Position', posE_main);
    axEcb   = axes(fig, 'Units', 'normalized', 'Position', posE_cb);

    plotPanelOmedaHeatmap_NoTaxa( ...
        axEmain, ...
        axEcb, ...
        O, ...
        omedaMethodNames);

    addPanelLabel(fig, posE_main, 'E');

    exportgraphics(fig, outBase + ".png", ...
        'Resolution', pngDPI, ...
        'BackgroundColor', 'white');

    exportgraphics(fig, outBase + ".pdf", ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

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

function plotPanelLoadingHeatmap(axMain, axTaxa, axCB, L_plot, taxaNames, colNames, SigMask)

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

    drawHeatmapGrid(axMain, nRows, nCols, 0.25);

    addHeatmapNumbersAbs(axMain, L_plot, 2, 0.55, 7.2);

    addStarMarkers(axMain, SigMask, 7.2, 0.23);

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
    cb.Label.String = 'C-ASCA loadings';
    cb.Label.FontName = 'Arial';
    cb.Label.FontSize = 9;
    cb.Label.FontWeight = 'normal';

    axis(axCB, 'off');
end

%% ============================================================
% PANEL E: oMEDA HEATMAP WITHOUT TAXA LABELS
%% ============================================================

function plotPanelOmedaHeatmap_NoTaxa(axMain, axCB, O, methodNames)

    methodNames = string(methodNames(:))';

    nRows = size(O,1);
    nCols = size(O,2);

    imagesc(axMain, O);
    clim(axMain, [-1 1]);
    colormap(axMain, blueWhiteRedCMap(256));

    axMain.XTick = 1:nCols;
    axMain.XTickLabel = cellstr(methodNames);
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

    drawHeatmapGrid(axMain, nRows, nCols, 0.25);

    addHeatmapNumbersAbs(axMain, O, 2, 0.65, 6.3);

    cb = colorbar(axMain);
    cb.Units = 'normalized';
    cb.Position = axCB.Position;

    cb.FontName = 'Arial';
    cb.FontSize = 8;
    cb.LineWidth = 0.5;
    cb.Box = 'off';
    cb.Ticks = [-1 0 1];
    cb.TickLabels = {'-1','0','1'};
    cb.Label.String = 'oMEDA over PCA';
    cb.Label.FontName = 'Arial';
    cb.Label.FontSize = 9;
    cb.Label.FontWeight = 'normal';

    axis(axCB, 'off');
end

%% ============================================================
% LOCAL FUNCTION:
% Draw heatmap grid
%% ============================================================

function drawHeatmapGrid(ax, nRows, nCols, lineWidth)

    if nargin < 4
        lineWidth = 0.25;
    end

    for r = 0.5:1:(nRows + 0.5)
        plot(ax, [0.5 nCols + 0.5], [r r], '-', ...
            'Color', [0.88 0.88 0.88], ...
            'LineWidth', lineWidth);
    end

    for c = 0.5:1:(nCols + 0.5)
        plot(ax, [c c], [0.5 nRows + 0.5], '-', ...
            'Color', [0.88 0.88 0.88], ...
            'LineWidth', lineWidth);
    end
end

%% ============================================================
% LOCAL FUNCTION:
% Add heatmap numbers as absolute values
%
% - Displays absolute values only
% - Converts 0.00 and -0.00 to 0
% - Keeps heatmap colours signed because the matrix itself is unchanged
%% ============================================================

function addHeatmapNumbersAbs(ax, M, nDec, whiteThreshold, fontSize)

    [nRows, nCols] = size(M);

    roundTol = 0.5 * 10^(-nDec);

    for i = 1:nRows
        for j = 1:nCols

            val = M(i,j);

            if isnan(val)
                labelStr = 'NA';
                valForColor = 0;
            else
                valAbs = abs(val);
                valRounded = round(valAbs, nDec);

                if valRounded < roundTol
                    labelStr = '0';
                    valForColor = 0;
                else
                    labelStr = sprintf(['%.' num2str(nDec) 'f'], valRounded);
                    valForColor = valAbs;
                end
            end

            if valForColor > whiteThreshold
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
% Add significance stars
%% ============================================================

function addStarMarkers(ax, SigMask, fontSize, yOffset)

    if isempty(SigMask)
        return;
    end

    for i = 1:size(SigMask,1)
        for j = 1:size(SigMask,2)
            if SigMask(i,j)
                text(ax, j, i + yOffset, '*', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontName', 'Arial', ...
                    'FontWeight', 'bold', ...
                    'FontSize', fontSize, ...
                    'Color', [0 0 0], ...
                    'Clipping', 'on');
            end
        end
    end
end

%% ============================================================
% LOCAL FUNCTION:
% Normalize vector
%% ============================================================

function v = normalizeVector(v)

    v = v(:);
    nrm = norm(v);

    if nrm < eps
        warning('A vector with near-zero norm was detected. Returning zeros.');
        v = zeros(size(v));
    else
        v = v / nrm;
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