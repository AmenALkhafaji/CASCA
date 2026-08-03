%% ============================================================
%  FULL SCHUBERT C-ASCA + oMEDA SCRIPT + FINAL A-E PANEL
%
%  Input:
%    diarrhea.mat containing a table named Datos
%
%  Expected structure of Datos:
%    taxa columns, Methods, Groups, observation
%
%  Schubert dataset:
%    N = 237 observations per method
%    Groups: H and CDI
%
%  Method order:
%     1 Raw Data
%     2 CLR+1
%     3 CLR-BMR
%     4 TSS
%     5 Rarefaction
%     6 CSS
%     7 edgeR-TMM
%     8 DESeq2
%     9 ALDEx2
%    10 ANCOM
%
%  Panel E contains the ten method outputs plus the
%  C-ASCA consensus output, giving 11 columns.
%% ============================================================
clear all ;
clearvars;
close all force;
clc;
drawnow;

java.lang.System.gc


%% ============================================================
% SETTINGS
%% ============================================================

analysisLevel = "schubert";

dataFileName = "diarrhea.mat";

outDir = fullfile( ...
    pwd, ...
    "Final_Schubert_CASCA_oMEDA_Figures" ...
);

if ~exist(outDir, "dir")
    mkdir(outDir);
end

exportResolution = 1200;

N = 237;
expectedMethods = 10;
expectedRowsPerMethod = N;
expectedTotalRows = expectedMethods * expectedRowsPerMethod;

methodOrder = [ ...
    "Raw Data", ...
    "CLR+1", ...
    "CLR-BMR", ...
    "TSS", ...
    "Rarefaction", ...
    "CSS", ...
    "edgeR-TMM", ...
    "DESeq2", ...
    "ALDEx2", ...
    "ANCOM" ...
];

groupOrder = ["Control", "Case"];

%% ============================================================
% LOAD SCHUBERT DATA
%% ============================================================

dataFile = load(dataFileName);

if ~isfield(dataFile, "Datos")
    error("%s does not contain a variable named Datos.", dataFileName);
end

Datos = dataFile.Datos;

if ~istable(Datos)
    error("Datos must be a MATLAB table.");
end

requiredMetadata = ["Methods", "Groups", "observation"];

if ~all(ismember(requiredMetadata, string(Datos.Properties.VariableNames)))
    error( ...
        "Datos must contain the metadata columns Methods, Groups, and observation." ...
    );
end

fprintf( ...
    "Loaded Schubert Datos dimensions: %d rows x %d columns\n", ...
    height(Datos), ...
    width(Datos) ...
);

if height(Datos) ~= expectedTotalRows
    error( ...
        ["Datos contains %d rows, but 10 methods with %d Schubert " ...
         "observations per method require %d rows."], ...
        height(Datos), ...
        expectedRowsPerMethod, ...
        expectedTotalRows ...
    );
end

%% ============================================================
% BUILD C-ASCA INPUT
%% ============================================================

metadataColumns = ismember( ...
    string(Datos.Properties.VariableNames), ...
    requiredMetadata ...
);

taxaTable = Datos(:, ~metadataColumns);
load diarrhea.mat;
taxaNames = string(lab);
nTaxa = numel(taxaNames);

X_raw_all_methods = table2array(taxaTable);

if ~isnumeric(X_raw_all_methods)
    error("All taxa columns in Datos must be numeric.");
end

if any(isnan(X_raw_all_methods), "all")
    error("The taxa matrix contains missing values.");
end

if any(isinf(X_raw_all_methods), "all")
    error("The taxa matrix contains infinite values.");
end

X2 = double(X_raw_all_methods);

%% ============================================================
% METHOD CODES AND LABELS
%% ============================================================

MethodCode = double(Datos.Methods);

if any(isnan(MethodCode))
    error("Datos.Methods contains missing or non-numeric values.");
end

if any(~ismember(MethodCode, 1:10))
    invalidCodes = unique(MethodCode(~ismember(MethodCode, 1:10)));
    error( ...
        "Datos.Methods contains invalid method codes: %s", ...
        mat2str(invalidCodes(:)') ...
    );
end

methodCounts = arrayfun( ...
    @(k) sum(MethodCode == k), ...
    1:10 ...
);

if any(methodCounts ~= expectedRowsPerMethod)
    error( ...
        ["Each method code must occur exactly %d times. " ...
         "Observed counts are: %s"], ...
        expectedRowsPerMethod, ...
        mat2str(methodCounts) ...
    );
end

Method = strings(size(MethodCode));

for k = 1:expectedMethods
    Method(MethodCode == k) = methodOrder(k);
end

if any(strlength(Method) == 0)
    error("At least one method code could not be mapped to a method name.");
end

%% ============================================================
% SCHUBERT GROUP CODING
%
% H   = -1 = Control
% CDI =  1 = Case
%% ============================================================

Groups = string(Datos.Groups);
Group = nan(size(Groups));

Group(strcmpi(Groups, "H")) = -1;
Group(strcmpi(Groups, "CDI")) = 1;

recognizedGroups = ...
    strcmpi(Groups, "H") | ...
    strcmpi(Groups, "CDI");

if any(~recognizedGroups)
    unknownGroups = unique(Groups(~recognizedGroups));
    error( ...
        "Unrecognized Schubert group labels: %s", ...
        strjoin(cellstr(unknownGroups), ", ") ...
    );
end

fprintf( ...
    "Schubert group counts per method block: H = %d, CDI = %d\n", ...
    sum(Groups(1:N) == "H"), ...
    sum(Groups(1:N) == "CDI") ...
);

%% ============================================================
% OBSERVATION IDENTIFIERS AND ROW ORDER
%% ============================================================

observation = double(Datos.observation);

if any(isnan(observation))
    error("Datos.observation contains missing or non-numeric values.");
end

expectedMethodSequence = repelem((1:10)', expectedRowsPerMethod);

if ~isequal(MethodCode(:), expectedMethodSequence)
    error( ...
        ["Rows must be ordered as %d rows for method 1, followed by %d " ...
         "rows for method 2, continuing through method 10."], ...
        N, ...
        N ...
    );
end

expectedObservationSequence = repmat((1:N)', expectedMethods, 1);

if ~isequal(observation(:), expectedObservationSequence)
    error( ...
        "Datos.observation must run from 1 to %d within every method block.", ...
        N ...
    );
end

groupReferenceOriginal = Groups(1:N);

for methodID = 1:expectedMethods

    idx = MethodCode == methodID;

    if ~isequal(Groups(idx), groupReferenceOriginal)
        error( ...
            "The H/CDI sample order differs in method block %d.", ...
            methodID ...
        );
    end

end

fprintf("Method counts: %s\n", mat2str(methodCounts));

%% ============================================================
% ORIGINAL SEQUENCING DEPTH
%
% Sequencing depth is calculated from the raw-data block only
% and then repeated for all ten method blocks. It must not be
% recalculated from transformed CLR or normalized matrices.
%% ============================================================

rawBlock = X2(1:N, :);
rawDepth = sum(rawBlock, 2);

if any(~isfinite(rawDepth)) || any(rawDepth <= 0)
    error("The raw-data block contains invalid or zero library sizes.");
end

depth = repmat(rawDepth, expectedMethods, 1);

%% ============================================================
% BLOCK-WISE SSQ NORMALIZATION
%% ============================================================

X = X2;
Y = X;
d = depth;

num_parts = expectedMethods;
row_size = N;

if size(X, 1) ~= num_parts * row_size
    error( ...
        "The taxa matrix must contain %d rows.", ...
        num_parts * row_size ...
    );
end

for i = 1:num_parts

    idx = (i - 1) * row_size + 1 : i * row_size;

    part = X(idx, :);
    part2 = preprocess2D(part, "Preprocessing", 1);

    sum_sq = sum(part2(:).^2);

    if ~isfinite(sum_sq)
        error("Non-finite SSQ detected in method block %d.", i);
    end

    if sum_sq > 0
        part = part / sqrt(sum_sq);
    end

    Y(idx, :) = part;
    d(idx) = tiedrank(rawDepth);

end

X = Y;
depth = d;

%% ============================================================
% RUN C-ASCA / PARGLM
%% ============================================================

F = [MethodCode, Group, observation, depth];

[T, cascaStruct] = parglm( ...
    X, ...
    F, ...
    "Preprocessing", 1, ...
    "Model", [1 2], ...
    "Ordinal", [0 0 0 1], ...
    "Random", [0 0 1 0], ...
    "Nested", [2 3] ...
);

disp("Schubert C-ASCA with SSQ block scaling completed.");

if height(T) >= 5
    T.Source(2:5) = { ...
        "Method", ...
        "Group", ...
        "observation", ...
        "depth" ...
    };
end

disp(T);

ascao = asca(cascaStruct);

%% ============================================================
% PARGLMMC SIGNIFICANCE FOR PANEL D
%% ============================================================

alpha = 0.01;
nPermMC = 10;

[T_mc, parglmMCStruct] = parglmMC( ...
    X, ...
    F, ...
    "Preprocessing", 1, ...
    "Model", [1 2], ...
    "Ordinal", [0 0 0 1], ...
    "Random", [0 0 1 0], ...
    "Nested", [2 3], ...
    "Permutations", nPermMC, ...
    "Mtc", 3 ...
);

disp("parglmMC table:");
disp(T_mc);

if ~isfield(parglmMCStruct, "p")
    error("parglmMCStruct does not contain field p.");
end

if size(parglmMCStruct.p, 1) ~= nTaxa
    error( ...
        "parglmMCStruct.p has %d rows, but %d taxa were detected.", ...
        size(parglmMCStruct.p, 1), ...
        nTaxa ...
    );
end

if size(parglmMCStruct.p, 2) < 2
    error( ...
        "parglmMCStruct.p must contain at least Method and Group columns." ...
    );
end

p_disagreement = parglmMCStruct.p(:, 1);
p_consensus = parglmMCStruct.p(:, 2);

SigLoad = [ ...
    p_disagreement(:) <= alpha, ...
    p_consensus(:) <= alpha ...
];

fprintf("\nparglmMC p-value diagnostics\n");
fprintf( ...
    "Disagreement: min %.3g | median %.3g | max %.3g | n <= %.3f: %d/%d\n", ...
    min(p_disagreement, [], "omitnan"), ...
    median(p_disagreement, "omitnan"), ...
    max(p_disagreement, [], "omitnan"), ...
    alpha, ...
    sum(p_disagreement <= alpha), ...
    numel(p_disagreement) ...
);
fprintf( ...
    "Consensus: min %.3g | median %.3g | max %.3g | n <= %.3f: %d/%d\n", ...
    min(p_consensus, [], "omitnan"), ...
    median(p_consensus, "omitnan"), ...
    max(p_consensus, [], "omitnan"), ...
    alpha, ...
    sum(p_consensus <= alpha), ...
    numel(p_consensus) ...
);

%% ============================================================
% PLOTTING LABELS
%% ============================================================

G = strings(size(Group));
G(Group == -1) = "Control";
G(Group == 1) = "Case";

Mtr = Method;

%% ============================================================
% SCORE PLOT A: METHOD-RELATED DISAGREEMENT
%% ============================================================

plotScores1D_BoxSquare_NoLegend( ...
    ascao.factors{1}.matrix, ...
    Mtr, ...
    1, ...
    methodOrder, ...
    fullfile(outDir, "Factor1_disagreement_scores_by_method"), ...
    exportResolution ...
);

%% ============================================================
% SCORE PLOT B: GROUP-RELATED CONSENSUS
%% ============================================================

plotScores1D_BoxSquare_NoLegend( ...
    ascao.factors{2}.matrix, ...
    G, ...
    1, ...
    groupOrder, ...
    fullfile(outDir, "Factor2_consensus_scores_by_group"), ...
    exportResolution ...
);

%% ============================================================
% SCORE PLOT C: COMBINED MODEL
%% ============================================================

if ~isfield(ascao, "interactions") || isempty(ascao.interactions)
    error( ...
        "No interaction matrix was found in ascao.interactions." ...
    );
end

Mcombined = ...
    ascao.factors{1}.matrix + ...
    ascao.factors{2}.matrix + ...
    ascao.interactions{1}.matrix;

plotCombinedModelScores_BoxSquare( ...
    Mcombined, ...
    Mtr, ...
    G, ...
    methodOrder, ...
    groupOrder, ...
    fullfile(outDir, "combined_model_scores"), ...
    exportResolution ...
);

%% ============================================================
% PANEL D: C-ASCA LOADING HEATMAP
%% ============================================================

a = ascao.factors{2}.scores(1);
b = Group(1);

if sign(a) == sign(b)
    load_consensus = ascao.factors{2}.loads(:, 1);
else
    load_consensus = -ascao.factors{2}.loads(:, 1);
end

load_disagreement = ascao.factors{1}.loads(:, 1);

L = [ ...
    load_disagreement(:), ...
    load_consensus(:) ...
];

meanCenterDisplayedValues = true;

if meanCenterDisplayedValues
    L_plot = L - mean(L, 1, "omitnan");
else
    L_plot = L;
end

tol = 5e-4;
L_plot(abs(L_plot) < tol) = 0;

%% ============================================================
% SELECT THE TOP 25 CONSENSUS TAXA
%
% Taxa are ranked by the absolute magnitude of the oriented,
% uncentred consensus loading. The same row indices are then used
% in Panels D and E, ensuring exact correspondence between them.
%% ============================================================

nTopConsensusTaxa = min(25, nTaxa);

[~, consensusRank] = sort( ...
    abs(load_consensus), ...
    "descend", ...
    "MissingPlacement", "last" ...
);

topConsensusIdx = consensusRank(1:nTopConsensusTaxa);

L_plot_top25 = L_plot(topConsensusIdx, :);
taxaNames_top25 = taxaNames(topConsensusIdx);
SigLoad_top25 = SigLoad(topConsensusIdx, :);

fprintf("\nTop %d taxa selected by absolute consensus loading:\n", ...
    nTopConsensusTaxa);

for rankID = 1:nTopConsensusTaxa
    taxonID = topConsensusIdx(rankID);
    fprintf( ...
        "%2d. %-35s consensus loading = % .6f\n", ...
        rankID, ...
        char(taxaNames(taxonID)), ...
        load_consensus(taxonID) ...
    );
end

plotTwoColumnLoadingHeatmap( ...
    L_plot_top25, ...
    taxaNames_top25, ...
    ["Disagreement", "Consensus"], ...
    SigLoad_top25, ...
    fullfile(outDir, "loadings_heatmap_top25_consensus"), ...
    exportResolution ...
);

%% ============================================================
% PANEL E: oMEDA OVER PCA
%
% Ten method outputs + one C-ASCA consensus output
%% ============================================================

rowsPerMethod = N;
numberOfMethods = expectedMethods;

X_oMEDA = X2;
groupReference = Group(1:N);

rj = cell(1, numberOfMethods);

for methodID = 1:numberOfMethods

    idx = (methodID - 1) * N + 1 : methodID * N;
    X_method = X_oMEDA(idx, :);

    rj{methodID} = omedaPca( ...
        X_method, ...
        1, ...
        X_method, ...
        groupReference, ...
        "Preprocessing", 1 ...
    );

end

%% ============================================================
% C-ASCA CONSENSUS DATA
%
% B + C(B), followed by averaging corresponding observations
% across the ten normalization-method blocks.
%% ============================================================

if numel(ascao.factors) < 3
    error( ...
        "The fitted C-ASCA object does not contain the nested observation factor." ...
    );
end

X_consensus_repeated = ...
    ascao.factors{2}.matrix + ...
    ascao.factors{3}.matrix;

nVariables = size(X_consensus_repeated, 2);

if size(X_consensus_repeated, 1) ~= N * numberOfMethods
    error( ...
        ["Expected the repeated consensus matrix to contain %d rows, " ...
         "but it contains %d rows."], ...
        N * numberOfMethods, ...
        size(X_consensus_repeated, 1) ...
    );
end

X_consensus_3D = reshape( ...
    X_consensus_repeated, ...
    N, ...
    numberOfMethods, ...
    nVariables ...
);

X_consensus = mean( ...
    X_consensus_3D, ...
    2, ...
    "omitnan" ...
);

X_consensus = reshape( ...
    X_consensus, ...
    N, ...
    nVariables ...
);

Group_consensus = Group(1:N);

rj_consensus = omedaPca( ...
    X_consensus, ...
    1, ...
    X_consensus, ...
    Group_consensus, ...
    "Preprocessing", 1 ...
);

%% ============================================================
% NORMALIZE oMEDA VECTORS AND BUILD PANEL E MATRIX
%% ============================================================

for methodID = 1:numberOfMethods
    rj{methodID} = normalizeVector(rj{methodID});
end

rj_consensus = normalizeVector(rj_consensus);

O = zeros(nTaxa, numberOfMethods + 1);

for methodID = 1:numberOfMethods
    O(:, methodID) = rj{methodID}(:);
end

O(:, numberOfMethods + 1) = rj_consensus(:);
O(abs(O) < tol) = 0;

if size(O, 1) ~= nTaxa
    error("The oMEDA matrix row count does not match the taxa names.");
end

omedaMethodNames = [ ...
    methodOrder, ...
    "C-ASCA" ...
];

%% ============================================================
% PANEL E SIGNIFICANCE MASK
%% ============================================================

nDecPanelE = 2;
zeroTolPanelE = 0.5 * 10^(-nDecPanelE);

Sig_OMEDA_raw = repmat( ...
    p_consensus(:) <= alpha, ...
    1, ...
    numberOfMethods + 1 ...
);

Sig_OMEDA = ...
    Sig_OMEDA_raw & ...
    abs(O) >= zeroTolPanelE;

% Use exactly the same taxa and row order selected for Panel D.
O_top25 = O(topConsensusIdx, :);
Sig_OMEDA_top25 = Sig_OMEDA(topConsensusIdx, :);

plotOmedaHeatmapCompact( ...
    O_top25, ...
    Sig_OMEDA_top25, ...
    omedaMethodNames, ...
    fullfile(outDir, "oMEDA_heatmap_top25_consensus_taxa"), ...
    exportResolution ...
);

%% ============================================================
% FINAL A-E PUBLICATION PANEL
%% ============================================================

createFinal_CASCA_oMEDA_Panel( ...
    ascao.factors{1}.matrix, ...
    ascao.factors{2}.matrix, ...
    Mcombined, ...
    Mtr, ...
    G, ...
    methodOrder, ...
    groupOrder, ...
    L_plot_top25, ...
    taxaNames_top25, ...
    SigLoad_top25, ...
    O_top25, ...
    Sig_OMEDA_top25, ...
    omedaMethodNames, ...
    fullfile(outDir, "Schubert_CASCA_oMEDA_final_panel"), ...
    analysisLevel ...
);

fprintf( ...
    "\nAll Schubert scores, heatmaps, and the final A-E panel were saved in:\n%s\n\n", ...
    outDir ...
);

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

%     exportgraphics(fig, outBase + ".png", ...
%         'Resolution', exportResolution, ...
%         'BackgroundColor', 'white');

%     exportgraphics(fig, outBase + ".pdf", ...
%         'ContentType', 'vector', ...
%         'BackgroundColor', 'white');

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

%     exportgraphics(figComb, outBase + ".png", ...
%         'Resolution', exportResolution, ...
%         'BackgroundColor', 'white');

%     exportgraphics(figComb, outBase + ".pdf", ...
%         'ContentType', 'vector', ...
%         'BackgroundColor', 'white');

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
    cellH = 0.32;

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

%     exportgraphics(figLoad, outBase + ".png", ...
%         'Resolution', exportResolution, ...
%         'BackgroundColor', 'white');

%     exportgraphics(figLoad, outBase + ".pdf", ...
%         'ContentType', 'vector', ...
%         'BackgroundColor', 'white');

    close(figLoad);
end

%% ============================================================
% LOCAL FUNCTION:
% Individual compact oMEDA heatmap WITHOUT taxa labels
%% ============================================================

function plotOmedaHeatmapCompact(O, SigMask, methodNames, outBase, exportResolution)

    methodNames = string(methodNames(:))';

    [nRows, nCols] = size(O);

    if nCols ~= numel(methodNames)
        error('Number of columns in O must match methodNames.');
    end

    if ~isequal(size(SigMask), size(O))
        error('SigMask must have the same size as O.');
    end

    showValues = true;
    nDec = 2;

    cellW = 0.48;
    cellH = 0.32;

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

    % Panel E stars intentionally removed.
    % addStarMarkersSuppressZero(axOmeda, SigMask, O, nDec, 8.0, 0.23);

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

%     exportgraphics(figOmeda, outBase + ".png", ...
%         'Resolution', exportResolution, ...
%         'BackgroundColor', 'white');
% % 
%     exportgraphics(figOmeda, outBase + ".pdf", ...
%         'ContentType', 'vector', ...
%         'BackgroundColor', 'white');

    close(figOmeda);
end

%% ============================================================
% FINAL PUBLICATION PANEL FUNCTION
% Fixed: A/B/C panel labels are no longer cropped at export.
%% ============================================================

function createFinal_CASCA_oMEDA_Panel( ...
    X_factor1, X_factor2, Mcombined, ...
    methodLabels, groupLabels, methodOrder, groupOrder, ...
    L_plot, taxaNames, SigLoad, O, SigOmeda, omedaMethodNames, outBase, analysisLevel)

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

    % ------------------------------------------------------------
    % Layout
    % Genus keeps the tall D/E heatmaps.
    % Phylum uses the shorter D/E heatmap height from the attached
    % phylum layout, because only five taxa otherwise create excessive
    % vertical stretching. Only the D/E heatmap row is changed for phylum.
    % ------------------------------------------------------------

    if lower(string(analysisLevel)) == "phylum"

        topY = 0.725;
        topH = 0.110;

        botY = 0.255;
        botH = 0.345;

    else

        topY = 0.835;
        topH = 0.115;

        botY = 0.090;
        botH = 0.680;

    end

    posA = [0.055  topY  0.245  topH];
    posB = [0.350  topY  0.210  topH];
    posC = [0.620  topY  0.325  topH];

    posD_taxa = [0.025  botY  0.135  botH];
    posD_main = [0.170  botY  0.180  botH];
    posD_cb   = [0.364  botY  0.018  botH];

    posE_main = [0.455  botY  0.420  botH];
    posE_cb   = [0.892  botY  0.018  botH];

    % ------------------------------------------------------------
    % Panel A: method-related disagreement scores
    % ------------------------------------------------------------

    axA = axes(fig, 'Units', 'normalized', 'Position', posA);
    plotPanelScore1D(axA, X_factor1, methodLabels, 1, methodOrder);
    addPanelLabel(fig, posA, 'A');

    % ------------------------------------------------------------
    % Panel B: group-related consensus scores
    % ------------------------------------------------------------

    axB = axes(fig, 'Units', 'normalized', 'Position', posB);
    plotPanelScore1D(axB, X_factor2, groupLabels, 1, groupOrder);
    addPanelLabel(fig, posB, 'B');

    % ------------------------------------------------------------
    % Panel C: combined model scores
    % ------------------------------------------------------------

    axC = axes(fig, 'Units', 'normalized', 'Position', posC);
    plotPanelCombinedScores(axC, Mcombined, methodLabels, groupLabels, methodOrder, groupOrder);
    addPanelLabel(fig, posC, 'C');

    % ------------------------------------------------------------
    % Panel D: C-ASCA loading heatmap
    % ------------------------------------------------------------

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

    addPanelLabel(fig, ...
        [posD_taxa(1), posD_main(2), posD_taxa(3) + posD_main(3), posD_main(4)], ...
        'D');

    % ------------------------------------------------------------
    % Panel E: oMEDA heatmap
    % ------------------------------------------------------------

    axEmain = axes(fig, 'Units', 'normalized', 'Position', posE_main);
    axEcb   = axes(fig, 'Units', 'normalized', 'Position', posE_cb);

    plotPanelOmedaHeatmap_NoTaxa(axEmain, axEcb, O, SigOmeda, omedaMethodNames);
    addPanelLabel(fig, posE_main, 'E');

    % ------------------------------------------------------------
    % Export
    % ------------------------------------------------------------

    drawnow;

%     exportgraphics(fig, outBase + ".png", ...
%         'Resolution', pngDPI, ...
%         'BackgroundColor', 'white');

    exportgraphics(fig, outBase + ".pdf", ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');

    close(fig);
end

%% ============================================================
% PANEL LABEL
% Fixed: labels are placed inside the panel boundary, not above it.
%% ============================================================

function addPanelLabel(fig, pos, txt)

    labelW = 0.045;
    labelH = 0.028;

    % Put the panel label at the true upper-left of each panel group.
    % The label is placed just above the panel and slightly to the left,
    % not inside the plotting area. This makes A/B/C/D/E read as panel
    % labels rather than data annotations.
    x = pos(1) - 0.030;
    y = pos(2) + pos(4) + 0.006;

    % Safety bounds to prevent cropping during PNG/PDF export.
    if x < 0.004
        x = 0.004;
    end

    if y + labelH > 0.996
        y = 0.996 - labelH;
    end

    annotation(fig, 'textbox', ...
        [x, y, labelW, labelH], ...
        'String', txt, ...
        'LineStyle', 'none', ...
        'FontName', 'Arial', ...
        'FontWeight', 'bold', ...
        'FontSize', 13, ...
        'Color', 'k', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
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
    xtickangle(ax, 45);

    ax.FontName = 'Arial';
    ax.FontSize = 5.5;
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

function plotPanelOmedaHeatmap_NoTaxa(axMain, axCB, O, SigMask, methodNames)

    methodNames = string(methodNames(:))';

    nRows = size(O,1);
    nCols = size(O,2);

    if ~isequal(size(SigMask), size(O))
        error('SigMask must have the same size as O.');
    end

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

    % Panel E stars intentionally removed.
    % addStarMarkersSuppressZero(axMain, SigMask, O, 2, 7.2, 0.23);

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

    [nRows, nCols] = size(SigMask);

    if isempty(SigMask)
        return;
    end

    if size(SigMask,1) ~= nRows
        % This guard is intentionally permissive because the axis data
        % dimensions are not directly available here. Dimension mismatch
        % will usually be caught before plotting.
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
% Add significance stars but suppress stars where the displayed value is 0
%
% Used for Panel E oMEDA heatmaps.
%% ============================================================

function addStarMarkersSuppressZero(ax, SigMask, ValueMatrix, nDec, fontSize, yOffset)

    if isempty(SigMask)
        return;
    end

    if ~isequal(size(SigMask), size(ValueMatrix))
        error('SigMask and ValueMatrix must have the same size.');
    end

    zeroTol = 0.5 * 10^(-nDec);

    [nRows, nCols] = size(SigMask);

    for i = 1:nRows
        for j = 1:nCols
            if SigMask(i,j) && abs(ValueMatrix(i,j)) >= zeroTol
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

        case "CSS"
            rgb = [0.3010 0.7450 0.9330];

        case "edgeR-TMM"
            rgb = [0.8500 0.3250 0.0980];

        case "DESeq2"
            rgb = [0.4660 0.6740 0.1880];

        case "ALDEx2"
            rgb = [0.6350 0.0780 0.1840];

        case "ANCOM"
            rgb = [0.2500 0.2500 0.2500];

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