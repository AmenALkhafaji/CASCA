%% ============================================================
%  General KW + BH-FDR + confusion matrix vs ground truth
%  Works for both PHYLUM and GENUS
%
%  Required variables:
%       X
%       ascao
%       GTP1   for phylum ground truth
%       GGT    for genus ground truth
%
%  Data structure:
%       X(1:300,:)      = Raw data
%       X(301:600,:)    = TSS
%       X(601:900,:)    = Rarefaction
%       X(901:1200,:)   = CLR
%
%  C-ASCA:
%       ascao.factors{2}.matrix
%       optionally ascao.factors{3}.matrix and ascao.factors{4}.matrix
%
%  Ground truth:
%       Phylum: GT = GTP1(:,2) - GTP1(:,1)
%       Genus : GT = GGT(2,:)  - GGT(1,:)
%% ============================================================

clear resultsTable kwTable taxonTable

clc;

%% ---------- Choose analysis level ----------
%analysisLevel = 'phylum';
 analysisLevel = 'genus';

%% ---------- Basic settings ----------
N = 300;
alpha = 0.05;
tolGT = 1e-12;

%% ---------- Select ground truth ----------
switch lower(string(analysisLevel))

    case "phylum"

        % Expected GTP1 format:
        % rows = phyla
        % column 1 = control
        % column 2 = case
        GT = GTP1(:,2) - GTP1(:,1);
        GT = GT(:)';   % force 1 x taxa

        tableSuffix = 'phylum';
        captionLevel = 'Phylum';
        expectedTaxa = 5;

    case "genus"

        % Expected GGT format:
        % row 1 = control
        % row 2 = case
        % columns = genera
        GT = GGT(:,2) - GGT(:,1);
        GT = GT(:)';   % force 1 x taxa

        tableSuffix = 'genus';
        captionLevel = 'Genus';
        expectedTaxa = 36;

    otherwise
        error('analysisLevel must be either ''phylum'' or ''genus''.');
end

nTaxa = numel(GT);

if nTaxa ~= expectedTaxa
    warning('%s analysis expected %d taxa, but ground truth has %d values.', ...
        captionLevel, expectedTaxa, nTaxa);
end

trueDA  = abs(GT) > tolGT;
trueDir = sign(GT);

trueDA  = trueDA(:)';     % force 1 x taxa
trueDir = trueDir(:)';    % force 1 x taxa

nTrueDA = sum(trueDA, 'all');

fprintf('\n====================================================\n');
fprintf('%s-level analysis\n', captionLevel);
fprintf('Number of true DA taxa = %d out of %d\n', nTrueDA, nTaxa);
fprintf('====================================================\n');

if nTrueDA == nTaxa
    warning(['All taxa are nonzero in the ground-truth vector. ', ...
             'Therefore FP = 0 for every method if confusion is based only on detected-vs-true DA taxa. ', ...
             'Precision will not be informative.']);
end

%% ---------- Safety checks for X ----------
if size(X,1) < 4*N
    error('X has %d rows, but at least %d rows are required for four method blocks.', ...
        size(X,1), 4*N);
end

if size(X,1) > 4*N
    warning(['X has %d rows, but this script expects the first %d rows only ', ...
             'for Raw data, TSS, Rarefaction, and CLR. Extra rows are ignored.'], ...
             size(X,1), 4*N);
end

if size(X,2) ~= nTaxa
    error(['X has %d taxa columns, but the selected ground-truth vector has %d values. ', ...
           'Check analysisLevel and the active X matrix.'], ...
           size(X,2), nTaxa);
end

%% ---------- Extract method matrices ----------
X_raw  = X(1:N, :);
X_tss  = X(N+1:2*N, :);
X_rare = X(2*N+1:3*N, :);
X_clr  = X(3*N+1:4*N, :);

%% ---------- Build C-ASCA representation ----------
xB = ascao.factors{2}.matrix;

if numel(ascao.factors) >= 4

    xC = ascao.factors{3}.matrix;
    xD = ascao.factors{4}.matrix;

    if ~isequal(size(xB), size(xC), size(xD))
        error('ascao.factors{2}, ascao.factors{3}, and ascao.factors{4} do not have the same size.');
    end

    Xmb_keep = xB + xC + xD;

else

    warning('ascao.factors{3} or ascao.factors{4} not available. Using factor 2 only for C-ASCA.');
    Xmb_keep = xB;

end

% Case 1: already N x taxa
if size(Xmb_keep,1) == N && size(Xmb_keep,2) == nTaxa

    X_CASCA = Xmb_keep;

% Case 2: stacked as N*K x taxa
elseif mod(size(Xmb_keep,1), N) == 0 && size(Xmb_keep,2) == nTaxa

    K = size(Xmb_keep,1) / N;
    X3 = reshape(Xmb_keep, N, K, nTaxa);
    X_CASCA = squeeze(mean(X3, 2, 'omitnan'));

% Case 3: transposed N x taxa
elseif size(Xmb_keep,2) == N && size(Xmb_keep,1) == nTaxa

    X_CASCA = Xmb_keep';

% Case 4: transposed stacked taxa x N*K
elseif mod(size(Xmb_keep,2), N) == 0 && size(Xmb_keep,1) == nTaxa

    K = size(Xmb_keep,2) / N;
    Xtemp = Xmb_keep';
    X3 = reshape(Xtemp, N, K, nTaxa);
    X_CASCA = squeeze(mean(X3, 2, 'omitnan'));

else

    error(['Cannot reshape C-ASCA matrix.\n', ...
           'size(Xmb_keep) = %d x %d\n', ...
           'Expected either N x taxa, N*K x taxa, taxa x N, or taxa x N*K.\n', ...
           'N = %d, nTaxa = %d'], ...
           size(Xmb_keep,1), size(Xmb_keep,2), N, nTaxa);

end

if size(X_CASCA,1) ~= N || size(X_CASCA,2) ~= nTaxa
    error('X_CASCA must be %d x %d, but current size is %d x %d.', ...
        N, nTaxa, size(X_CASCA,1), size(X_CASCA,2));
end

%% ---------- Methods ----------
methods = {'Raw data', 'TSS', 'Rarefaction', 'CLR', 'C-ASCA'};
dataCell = {X_raw, X_tss, X_rare, X_clr, X_CASCA};

nMethods = numel(methods);

%% ---------- Group labels ----------
% Assumption:
%   rows 1:150   = control
%   rows 151:300 = case
%
% Important:
%   This uses 0/1 coding only for the KW test and median direction.
%   effect = median(case) - median(control)

group = [zeros(150,1); ones(150,1)];

if numel(group) ~= N
    error('Group vector must have %d rows.', N);
end

%% ---------- Safety checks ----------
for m = 1:nMethods

    Xm = dataCell{m};

    if size(Xm,1) ~= N
        error('%s has %d rows. Expected %d rows.', methods{m}, size(Xm,1), N);
    end

    if size(Xm,2) ~= nTaxa
        error('%s has %d taxa columns, but ground truth has %d values.', ...
            methods{m}, size(Xm,2), nTaxa);
    end

end

%% ---------- Preallocate ----------
TP = zeros(nMethods,1);
FP = zeros(nMethods,1);
FN = zeros(nMethods,1);

Precision = zeros(nMethods,1);
Recall    = zeros(nMethods,1);
F1        = zeros(nMethods,1);

nSig = zeros(nMethods,1);
wrongDirection = zeros(nMethods,1);

allP = nan(nMethods, nTaxa);
allQ = nan(nMethods, nTaxa);
allSig = false(nMethods, nTaxa);
allEffect = nan(nMethods, nTaxa);
allPredDir = nan(nMethods, nTaxa);

%% ---------- Main analysis ----------
for m = 1:nMethods

    Xm = dataCell{m};

    pvals = nan(1,nTaxa);
    effect = nan(1,nTaxa);

    for j = 1:nTaxa

        xj = Xm(:,j);

        valid = ~isnan(xj) & ~isnan(group);
        xj_valid = xj(valid);
        group_valid = group(valid);

        if numel(unique(group_valid)) < 2
            pvals(j) = NaN;
            effect(j) = NaN;
            continue;
        end

        if numel(unique(xj_valid)) <= 1
            pvals(j) = 1;
            effect(j) = 0;
            continue;
        end

        % Kruskal-Wallis test: control vs case
        pvals(j) = kruskalwallis(xj_valid, group_valid, 'off');

        % Estimated direction: median case - median control
        medControl = median(xj_valid(group_valid == 0), 'omitnan');
        medCase    = median(xj_valid(group_valid == 1), 'omitnan');

        effect(j) = medCase - medControl;

    end

    % BH-FDR correction
    qvals = bh_fdr_local(pvals);

    % Significant taxa
    sig = qvals <= alpha;

    % Predicted direction
    predDir = sign(effect);

    % Force all vectors to same orientation
    sig     = sig(:)';
    predDir = predDir(:)';
    trueDA  = trueDA(:)';
    trueDir = trueDir(:)';

    if numel(sig) ~= nTaxa || numel(trueDA) ~= nTaxa || ...
       numel(predDir) ~= nTaxa || numel(trueDir) ~= nTaxa
        error('Vector length mismatch before confusion calculation.');
    end

    % Directional agreement with ground truth
    dirOK = predDir == trueDir;

    % Confusion matrix
    %
    % TP = significant and truly DA
    % FP = significant but not truly DA
    % FN = not significant but truly DA
    %
    % Wrong direction is reported separately.
    TP(m) = sum(sig & trueDA, 'all');
    FP(m) = sum(sig & ~trueDA, 'all');
    FN(m) = sum(~sig & trueDA, 'all');

    if TP(m) + FP(m) > 0
        Precision(m) = TP(m) / (TP(m) + FP(m));
    else
        Precision(m) = NaN;
    end

    if TP(m) + FN(m) > 0
        Recall(m) = TP(m) / (TP(m) + FN(m));
    else
        Recall(m) = NaN;
    end

    if Precision(m) + Recall(m) > 0
        F1(m) = 2 * Precision(m) * Recall(m) / (Precision(m) + Recall(m));
    else
        F1(m) = NaN;
    end

    nSig(m) = sum(sig, 'all');

    % Significant true taxa but opposite estimated direction
    wrongDirection(m) = sum(sig & trueDA & ~dirOK, 'all');

    allP(m,:) = pvals(:)';
    allQ(m,:) = qvals(:)';
    allSig(m,:) = sig(:)';
    allEffect(m,:) = effect(:)';
    allPredDir(m,:) = predDir(:)';

end

%% ---------- Results table ----------
resultsTable = table( ...
    methods(:), TP, FP, FN, Precision, Recall, F1, wrongDirection, ...
    'VariableNames', {'Method','TP','FP','FN','Precision','Recall','F1','WrongDirection'} ...
);

disp(' ');
disp('Confusion-matrix results against ground truth');
disp(resultsTable);

%% ---------- KW discovery count table ----------
sigText = strings(nMethods,1);

for m = 1:nMethods
    sigText(m) = sprintf('%d/%d', nSig(m), nTaxa);
end

kwTable = table(methods(:), sigText, ...
    'VariableNames', {'Method','q_le_0_05'} ...
);

disp(' ');
disp('Kruskal-Wallis BH-FDR discovery counts');
disp(kwTable);

%% ---------- Taxon-level diagnostic table ----------
taxonTable = table;
taxonTable.TaxonIndex = (1:nTaxa)';
taxonTable.GroundTruth = GT(:);
taxonTable.TrueDA = trueDA(:);
taxonTable.TrueDirection = trueDir(:);

for m = 1:nMethods

    cleanName = matlab.lang.makeValidName(methods{m});

    taxonTable.([cleanName '_p'])       = allP(m,:)';
    taxonTable.([cleanName '_q'])       = allQ(m,:)';
    taxonTable.([cleanName '_sig'])     = allSig(m,:)';
    taxonTable.([cleanName '_effect'])  = allEffect(m,:)';
    taxonTable.([cleanName '_predDir']) = allPredDir(m,:)';

end

disp(' ');
disp('Taxon-level diagnostic table');
disp(taxonTable);

%% ---------- Save tables ----------
outDir = fullfile(pwd, ['KW_confusion_' tableSuffix]);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

writetable(resultsTable, fullfile(outDir, ['confusion_' tableSuffix '.csv']));
writetable(kwTable,      fullfile(outDir, ['kw_counts_' tableSuffix '.csv']));
writetable(taxonTable,   fullfile(outDir, ['taxon_diagnostics_' tableSuffix '.csv']));

fprintf('\nSaved results in folder:\n%s\n', outDir);

%% ---------- Print LaTeX confusion table ----------
fprintf('\n\n%% ============================================================\n');
fprintf('%% LaTeX confusion-matrix table\n');
fprintf('%% ============================================================\n');

fprintf('\\begin{table}[H]\n');
fprintf('\\centering\n');
fprintf('\\caption{Confusion-matrix results (TP, FP, FN) and derived metrics (precision, recall, F1) with respect to the ground truth for the %s simulation. Wrong-direction detections indicate taxa detected as significant but with an estimated effect direction opposite to the simulated ground truth.}\n', lower(captionLevel));
fprintf('\\label{tab:%s_confusion}\n', tableSuffix);
fprintf('\\scriptsize\n');
fprintf('\\setlength{\\tabcolsep}{2pt}\n');
fprintf('\\renewcommand{\\arraystretch}{1.0}\n');
fprintf('\\begin{tabularx}{\\columnwidth}{@{}l *{7}{>{\\centering\\arraybackslash}X}@{}}\n');
fprintf('\\toprule\n');
fprintf('\\textbf{Method} & \\textbf{TP} & \\textbf{FP} & \\textbf{FN} & \\textbf{Precision} & \\textbf{Recall} & \\textbf{F1} & \\textbf{Wrong dir.} \\\\\n');
fprintf('\\midrule\n');

for m = 1:nMethods
    fprintf('%s & %d & %d & %d & %.2f & %.2f & %.2f & %d \\\\\n', ...
        methods{m}, TP(m), FP(m), FN(m), Precision(m), Recall(m), F1(m), wrongDirection(m));
end

fprintf('\\bottomrule\n');
fprintf('\\end{tabularx}\n');
fprintf('\\end{table}\n');

%% ---------- Print LaTeX KW table ----------
fprintf('\n\n%% ============================================================\n');
fprintf('%% LaTeX KW discovery-count table\n');
fprintf('%% ============================================================\n');

fprintf('\\begin{table}[H]\n');
fprintf('\\centering\n');
fprintf('\\caption{%s-level Kruskal--Wallis discovery counts after BH--FDR control. For each preprocessing output, a Kruskal--Wallis test was applied for the control versus case contrast and the resulting $p$-values were adjusted using the Benjamini--Hochberg procedure. Entries report the number of significant taxa out of the total (%d) at $q\\le 0.05$.}\n', captionLevel, nTaxa);
fprintf('\\label{tab:kw%s}\n', tableSuffix);
fprintf('\\scriptsize\n');
fprintf('\\renewcommand{\\arraystretch}{1.0}\n');
fprintf('\\setlength{\\tabcolsep}{4pt}\n');
fprintf('\\begin{tabularx}{\\columnwidth}{@{}X >{\\raggedleft\\arraybackslash}p{0.24\\columnwidth}@{}}\n');
fprintf('\\toprule\n');
fprintf('\\textbf{Method} & \\textbf{$q\\le 0.05$} \\\\\n');
fprintf('\\midrule\n');

for m = 1:nMethods
    fprintf('%s & %d/%d \\\\\n', methods{m}, nSig(m), nTaxa);
end

fprintf('\\bottomrule\n');
fprintf('\\end{tabularx}\n');
fprintf('\\end{table}\n');

%% ============================================================
%  Local function: Benjamini-Hochberg FDR
%% ============================================================
function q = bh_fdr_local(p)

    p = p(:);
    q = nan(size(p));

    valid = ~isnan(p);
    pv = p(valid);

    if isempty(pv)
        q = q(:)';
        return;
    end

    [pv_sorted, order] = sort(pv, 'ascend');
    m = numel(pv_sorted);

    q_sorted = pv_sorted .* m ./ (1:m)';

    % Enforce monotonicity
    for i = m-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end

    q_sorted(q_sorted > 1) = 1;

    q_valid = nan(size(pv));
    q_valid(order) = q_sorted;

    q(valid) = q_valid;

    q = q(:)';   % return row vector

end


eff_CASCA_check = nan(1,nTaxa);

for j = 1:nTaxa
    xj = X_CASCA(:,j);

    medControl = median(xj(group == 0), 'omitnan');
    medCase    = median(xj(group == 1), 'omitnan');

    eff_CASCA_check(j) = medCase - medControl;
end

GT_row = GT(:)';

r_casca_gt = corr(eff_CASCA_check(:), GT_row(:), 'Rows', 'complete');

fprintf('\nCorrelation between C-ASCA effect and ground truth = %.4f\n', r_casca_gt);

disp(table((1:nTaxa)', GT_row(:), eff_CASCA_check(:), ...
    sign(GT_row(:)), sign(eff_CASCA_check(:)), ...
    'VariableNames', {'Taxon','GT','CASCA_Effect','GT_Dir','CASCA_Dir'}));