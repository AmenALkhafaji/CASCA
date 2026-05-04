% X: your data matrix (samples x variables)
% group: a vector with group labels (samples x 1)
clc;
clear all;
close all;
 load genus.mat;
load GenusASCA.mat;
  %load phylum.mat;
% load PhylumASCA.mat;
X =table2array(Datos(:,1:(size(Datos,2)-3)));
depth=(sum(X,2));
Method = string(Datos.Methods);
Method(Method=='raw data')=1;
Method(Method=='TSS')=2;
Method(Method=='rarefaction')=3;
Method(Method=='CLR')=4;
 Method(Method=='rCLR')=4;
Method=double(Method);
Groups = Datos.Groups;  
%Groups = Datos.tags; % categorical
Group = ones(size(Groups));     % double
Group(Groups=="Healthy") = -1;
 Group(Groups=="sick")    =  1;

%  Group(Groups==1) = -1;
%  Group(Groups==2)    =  1;
observation =(Datos.observation);
Y=X;
num_parts=4;
row_size=size(X, 1)/4;
d=depth;

for i = 1:num_parts

    part = X((i-1)*row_size + 1:i*row_size, :);
    part2 = preprocess2D(part,'Preprocessing',1);

    % Compute the sum of squares
    sum_sq = sum(part2(:).^2);
    
    % Normalize the part
    if sum_sq > 0
        part = part / sqrt(sum_sq);
    end
Y((i-1)*row_size + 1:i*row_size, :)=part;
 d((i-1)*row_size + 1:i*row_size, :)=tiedrank(depth((i-1)*row_size + 1:i*row_size, :));

end
depth=(d);

X=Y;
F=[Method, Group, observation, depth];
[T, struct] = parglm(X, F, 'Preprocessing',1,'Model', [1 2],'Ordinal',[0 0 0 1],'Random',[0 0 1 0],'Nested',[2,3]);
disp("ASCA with SSQ")
T.Source(2:5) = { 'Method'  'Group' 'observation' 'depth'}
ascao = asca(struct);

model=ascao;
% 
G=string(Group);
G(G=="-1")="CTRL";
G(G=="1")="Case";


Mtr = string(Datos.Methods);
Mtr(Mtr=='raw data')="Raw Data";
Mtr(Mtr=='TSS')="TSS";
Mtr(Mtr=='rarefaction')="Rarefaction";
Mtr(Mtr=='CLR')="CLR";

% plotVec(sum(ascao.residuals.^2,2),'ObsClass', Mtr,'XYLabel',{'','Sq. Residuals'});
% plotVec(sum(ascao.residuals.^2,2),'ObsClass',G,'XYLabel',{'','Sq. Residuals'}); 
% plotVec(sum(ascao.residuals.^2,2),'ObsClass',depth,'XYLabel',{'','Sq. Residuals'}); 
% 
% %% Visualization
% 
% 
% %% factor = 1
% z2(ascao.factors{2}.loads>0)="Case";
% z2(ascao.factors{2}.loads<0)="Control";
% 
%   scoresPca(ascao.factors{1}.matrix,'PCs',1,'ObsTest',ascao.factors{1}.matrix,'Preprocessing',1,'PlotCal',false,'ObsClass',Mtr);
%  loadingsPca(ascao.factors{1},'PCs',1:2,'Preprocessing',1,'VarsClass',z2,'VarsLabel',lab);
% 
% 
% 
% 
% 
% 
% %% factor = 2
% z2(ascao.factors{2}.loads>0)="Case";
% z2(ascao.factors{2}.loads<0)="Control";
%  scoresPca(ascao.factors{2}.matrix,'PCs',1,'ObsTest',X,'Preprocessing',0,'PlotCal',false,'ObsClass',G);
%  loadingsPca(ascao.factors{2},'PCs',1,'Preprocessing',1,'VarsClass',z2,'VarsLabel',lab);
% 
% 
%  %% interaction
%  Ff=string(Datos.Methods);
%  M = ascao.factors{1}.matrix + ascao.factors{2}.matrix + ascao.interactions{1}.matrix;
% codeLevels = {};
% for i=1:size(F,1), codeLevels{i} = sprintf(':%s,:%s',Mtr(i),G(i));end;
% scoresPca(M,'PCs',1:2,'ObsTest',M,'Preprocessing',1,'PlotCal',false,'ObsClass',codeLevels);
% loadingsPca(M,'PCs',1,'Preprocessing',1,'VarsClass',lab,'VarsLabel',lab);
% legend(unique(codeLevels),'Location', 'eastoutside')
% 
% %% Error between oMEDA and ground truth
% % rj=omedaPca(ascao.factors{2}.matrix,[1], ascao.factors{2}.matrix, Group,'Preprocessing', 1);
% % GTP1=GGT;
% % GTP1 = GTP1 - mean(GTP1(:,1:2), 2) * ones(1, 3);
% % 
% % ground_truth = sign(GTP1(:,2)) .* (GTP1(:,2).^2)- sign(GTP1(:,1)) .* (GTP1(:,1).^2);
% % tot2 = rj/ norm(rj);
% % gt2 = ground_truth/ norm(ground_truth);
% % e2a = sum((tot2 - gt2).^2);
% % e2b = sum((tot2 + gt2).^2);
% % e2=min(e2a,e2b);
% % fprintf('Total Error  pca_omeda:%.4f\n', e2);
% 
% 
% 
% %% Error between loads and ground truth
% % rj=ascao.factors{2}.loads;
% % GTP1=GGT;
% %  GTP1 = GTP1 - mean(GTP1(:,1:2), 2) * ones(1, 3);
% % ground_truth = GTP1(:,2)-GTP1(:,1);
% % tot2 = rj/ norm(rj);
% % gt2 = ground_truth / norm(ground_truth);
% % e2a = sum((tot2 - gt2).^2);
% % e2b = sum((tot2 + gt2).^2);
% % e2=min(e2a,e2b);
% % fprintf('Total Error  load:%.4f\n', e2);
% 
% %% Consensus Factor export 
% 
% % N = 300;
% % 
% % xA  = ascao.factors{1}.matrix;
% % xAB = ascao.interactions{1}.matrix;
% % xB  = ascao.factors{2}.matrix;
% % xC  = ascao.factors{3}.matrix;    
% % xD  = ascao.factors{4}.matrix;
% % % Build a "consensus signal" (choose what you want to keep)
% % Xmb_keep = xA+xAB+xB + xC + xD;               % example: keep B + C(B) + D
% %                   K = size(Xmb_keep,1)/N;
% % X3 = reshape(Xmb_keep, N, K, []);
% 
% 
% % X_keep = squeeze(mean(X3,2,'omitnan'));   % N x p
% 
% 
% % Compute oMEDA vectors
% rj1 = omedaPca(X(1:300,:),    [1], X(1:300,:),    Group(1:300), 'Preprocessing', 1); % RawData
% rj2 = omedaPca(X(301:600,:),  [1], X(301:600,:),  Group(1:300), 'Preprocessing', 1); % TSS
% rj3 = omedaPca(X(601:900,:),  [1], X(601:900,:),  Group(1:300), 'Preprocessing', 1); % Rarefaction
% rj4 = omedaPca(X(901:1200,:), [1], X(901:1200,:), Group(1:300), 'Preprocessing', 1); % CLR
% rj5 = omedaPca(ascao.factors{2}.matrix, [1], ascao.factors{2}.matrix, Group, 'Preprocessing', 1); % C-ASCA
% 
% % Normalize safely to unit norm
% rj1 = normalizeVector(rj1);
% rj2 = normalizeVector(rj2);
% rj3 = normalizeVector(rj3);
% rj4 = normalizeVector(rj4);
% rj5 = normalizeVector(rj5);
% 
% 
% % Labels and method names
% lab = [
%     "Bifidobacterium"
%     "Butyricimonas"
%     "Odoribacter"
%     "Paraprevotella"
%     "Bacteroides"
%     "Parabacteroides"
%     "Prevotella"
%     "unknown1"
%     "unknown2"
%     "unknown3"
%     "Clostridium"
%     "unknown4"
%     "unknown5"
%     "Ruminococcus"
%     "Anaerostipes"
%     "Blautia"
%     "Coprococcus"
%     "Dorea"
%     "Lachnospira"
%     "Roseburia"
%     "unknown6"
%     "unknown7"
%     "unknown8"
%     "unknown9"
%     "Faecalibacterium"
%     "Oscillospira"
%     "Ruminococcus.1"
%     "Dialister"
%     "unknown10"
%     "unknown11"
%     "unknown12"
%     "Sutterella"
%     "Escherichia"
%     "Klebsiella"
%     "Haemophilus"
%     "Akkermansia"
% ];
% 
% taxaNames   = lab;
% methodNames = {'RawData', 'TSS', 'Rarefaction', 'CLR', 'C-ASCA'};
% 
% % Combine vectors
% dataMatrix = [rj1(:), rj2(:), rj3(:), rj4(:), rj5(:)];
% 
% % Sanity checks
% if size(dataMatrix,1) ~= numel(taxaNames)
%     error('Number of taxa labels (%d) does not match number of rows in dataMatrix (%d).', ...
%         numel(taxaNames), size(dataMatrix,1));
% end
% 
% if any(~isfinite(dataMatrix(:)))
%     error('dataMatrix contains NaN or Inf values. Check oMEDA outputs and normalization.');
% end
% 
% % Remove tiny numerical values to avoid showing -0.000
% tol = 5e-4;   % matched to 3 decimal places
% dataMatrix(abs(dataMatrix) < tol) = 0;
% 
% % Make taxa labels italic for display
% taxaNamesItalic = cellstr(compose('\\it{%s}', string(taxaNames)));
% 
% % Heatmap
% figure('Color', 'w', 'Position', [100 100 600 500]);
% h = heatmap(methodNames, taxaNamesItalic, dataMatrix, ...
%     'Colormap', redbluecmap());
% 
% % Display values rounded to 3 decimal places
% h.CellLabelFormat = '%.3f';
% 
% % Symmetric color limits around zero
% maxVal = max(abs(dataMatrix(:)));
% if maxVal > 0
%     h.ColorLimits = [-maxVal, maxVal];
% else
%     h.ColorLimits = [-1, 1];
% end
% 
% 
% % -------- Local helper functions --------
% function v = normalizeVector(v)
%     v = v(:);
%     nrm = norm(v);
%     if nrm < eps
%         warning('A vector with near-zero norm was detected. Returning zeros instead of dividing by zero.');
%         v = zeros(size(v));
%     else
%         v = v / nrm;
%     end
% end
% 
% function cmap = redbluecmap()
%     % Blue -> White -> Red
%     c1 = [0 0 1];
%     c2 = [1 1 1];
%     c3 = [1 0 0];
%     n = 64;
% 
%     cmap = [interp1([1 32 64], [c1(1) c2(1) c3(1)], 1:n)', ...
%             interp1([1 32 64], [c1(2) c2(2) c3(2)], 1:n)', ...
%             interp1([1 32 64], [c1(3) c2(3) c3(3)], 1:n)'];
% end
% % rj=omedaPca(X(1:300,:),[1],X(1:300,:),Group(1:300),'Preprocessing', 1)
