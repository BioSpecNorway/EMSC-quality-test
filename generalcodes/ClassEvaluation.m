function [YValH,Results] = ClassEvaluation(ClasModel,X,Levels,PLSPar,PlotIt)
% function [YHat] = ClassEvaluation(ClasModel,X,Levels,Vote,PlotIt)
% validation of classification model in ClasModel, plots and returns 
% classification results
% inputs:   ClasModel   hierarhical tree
%           X           saisir structure 
%           Levels      structure containing information about the tree
%                       (Names, NameRange, NumLevels)
%           PLSPar      parameters of PLSDA
%           PlotIt      two parameters 
%           PlotIt.Models for plotting confusion matrices; 
%           Plot.Levels for plotting confusion matrices for each level           
% outputs:  YValH       1xNum.Levels cells of predictions
%           Results     1xNum.Levels cells of prediction statistics
%
%% V.Tafintseva 10.05.2016
% 

YValH = cell(1,Levels.NumLevels);
which = (1:size(X.d,1))';
level = 1;

[YValH] = Classification(level,ClasModel,X,which,YValH,Levels,PLSPar);
YVal = X.i;

%% Validation for each node; plot results
if PlotIt.Models == 1
    PlotConfM4Validation(level,ClasModel,which,Levels,YVal,cell2mat(YValH));
end

%% Validation for all levels
Results = cell(1,Levels.NumLevels);
for i=1:Levels.NumLevels
   res = ClassResult(YVal(:,Levels.NameRange{i}),...
        YValH{i}(:,1:length(Levels.NameRange{i})));
    Results{i} = res{1};
end
% Plot results
if PlotIt.Levels == 1
    for i=1:Levels.NumLevels
        plotConfusionMat(Results{i},'val')
    end
end

end