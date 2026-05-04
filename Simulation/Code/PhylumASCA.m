

%%  PHYLUM WITHOUT ssbLOCK
clc;
clear all;
close all;
load PhylumASCA.mat;
X =table2array(Datos(:,1:(size(Datos,2)-3)));
Method =string(table2array(Datos(:,(size(Datos,2)-1)))); % method 
Group =string(table2array(Datos(:,(size(Datos,2)-2)))); % Group 
observation= string(table2array(Datos(:,(size(Datos,2))))); % observation
data=X;
Depth=sum(data,2);
F = [ Method Group  Depth observation];
[T, struct] = parglm(X, F(:,1:4), 'Preprocessing',0,'Permutations',1,'Model', [1 2],'Nested',[2 3]);  
T.Source(2:5) = { 'Method','Group ' 'Depth' 'observation'}
ascao = asca(struct);

scoresPca(ascao.factors{1}.matrix,'PCs',1:2,'ObsTest',X,'Preprocessing',0,'PlotCal',false,'ObsClass',Method);
legend(unique(Method), 'Location', 'eastoutside')
title('Scores of Method the data without blocks ');
scoresPca(microbiome_data,'PCs',1:2,"ObsTest",microbiome_data,'ObsClass',s,'ObsLabel',myemptylabel,'PlotCal', false);



title('Scores the data without blocks between group & methods');M = ascao.factors{1}.matrix + ascao.factors{2}.matrix + ascao.interactions{1}.matrix;
codeLevels = {};
for i=1:size(F,1), codeLevels{i} = sprintf(':%s,:%s',F(i,1),F(i,2));end;
scoresPca(M,'PCs',1:2,'ObsTest',X,'Preprocessing',0,'PlotCal',false,'ObsClass',codeLevels);
legend(unique(codeLevels),'Location', 'eastoutside')
plotVec(sum(ascao.residuals.^2,2),'ObsClass', Method,'XYLabel',{'','Sq. Residuals'});
plotVec(sum(ascao.residuals.^2,2),'ObsClass',Group,'XYLabel',{'','Sq. Residuals'}); 
plotVec(sum(ascao.residuals.^2,2),'ObsClass',Depth,'XYLabel',{'','Sq. Residuals'});

 %% PHYLUM WITH SSBLOCK 
Y=X;
num_parts=4;
row_size=size(X, 1)/4;
for i = 1:num_parts
    part = X((i-1)*row_size + 1:i*row_size, :);
    
    % Compute the sum of squares
    sum_sq = sum(part(:).^2);
    
    % Normalize the part
    if sum_sq > 0
        part = part / sqrt(sum_sq);
    end
Y((i-1)*row_size + 1:i*row_size, :)=part;
end

X=Y;
[T, struct] = parglm(X, F(:,1:4), 'Preprocessing',0,'Permutations',1,'Model', [1 2],'Nested',[2 3]);  
T.Source(2:5) = { 'Method','Group ' 'Depth' 'observation'}
ascao = asca(struct);

scoresPca(ascao.factors{1}.matrix,'PCs',1:2,'ObsTest',X,'Preprocessing',1,'PlotCal',false,'ObsClass',Method);
legend(unique(Method), 'Location', 'eastoutside')
title('Scores of Method the data with SSQ blocks ');


M = ascao.factors{1}.matrix + ascao.factors{2}.matrix + ascao.interactions{1}.matrix;
codeLevels = {};
for i=1:size(F,1), codeLevels{i} = sprintf(':%s,:%s',F(i,1),F(i,2));end;
scoresPca(M,'PCs',1:2,'ObsTest',X,'Preprocessing',0,'PlotCal',false,'ObsClass',codeLevels);
legend(unique(codeLevels),'Location', 'eastoutside')
title('Scores the data with SSQ blocks between group & methods');
plotVec(sum(ascao.residuals.^2,2),'ObsClass', Method,'XYLabel',{'','Sq. Residuals'});
plotVec(sum(ascao.residuals.^2,2),'ObsClass',Group,'XYLabel',{'','Sq. Residuals'}); 
plotVec(sum(ascao.residuals.^2,2),'ObsClass',Depth,'XYLabel',{'','Sq. Residuals'});

 %% PHYLUM WITH RANKTRANSFORMATION FOR ALL BLOCKS 
load PhylumASCA.mat;
X =table2array(Datos(:,1:(size(Datos,2)-3)));
X = rankTransform(X);
[T, struct] = parglm(X, F(:,1:4), 'Preprocessing',0,'Permutations',1,'Model', [1 2],'Nested',[2 3]);  
T.Source(2:5) = { 'Method','Group ' 'Depth' 'observation'}
ascao = asca(struct);

scoresPca(ascao.factors{1}.matrix,'PCs',1:2,'ObsTest',X,'Preprocessing',0,'PlotCal',false,'ObsClass',Method);
legend(unique(Method), 'Location', 'eastoutside')
title('Scores of Method the data RANKTRANSFORMATION  blocks ');
M = ascao.factors{1}.matrix + ascao.factors{2}.matrix + ascao.interactions{1}.matrix;
codeLevels = {};
for i=1:size(F,1), codeLevels{i} = sprintf(':%s,:%s',F(i,1),F(i,2));end;
scoresPca(M,'PCs',1:2,'ObsTest',X,'Preprocessing',0,'PlotCal',false,'ObsClass',codeLevels);
legend(unique(codeLevels),'Location', 'eastoutside')
title('Scores the data RANKTRANSFORMATION  blocks between group & methods');
plotVec(sum(ascao.residuals.^2,2),'ObsClass', Method,'XYLabel',{'','Sq. Residuals'});
plotVec(sum(ascao.residuals.^2,2),'ObsClass',Group,'XYLabel',{'','Sq. Residuals'}); 
plotVec(sum(ascao.residuals.^2,2),'ObsClass',Depth,'XYLabel',{'','Sq. Residuals'});


