function [YHat] = Classification(level,ClasModel,X,which,YHat,Levels,PLSPar)
% function [YHat] = Classification(level,ClasModel,X,which,YHat,Levels)
% iterative function, classifies samples from X according to a hierarhical 
% tree in ClasModel
% inputs:   level       number of the current level
%           ClasModel   hierarhical tree
%           X           saisir structure 
%           which       sample index on the current level
%           YHat        iteratively filled predictions
%           Levels      structure containing information about the tree
%           
% outputs:  YHat        predictions
%
%% V.Tafintseva 10.05.2016
%
if ~isempty(ClasModel.Model)
    Mod = ClasModel.Model;
    NX = size(X.d,1);
    YHatExact = [ones(NX,1),X.d]*[Mod.B0.d;Mod.B.d];
    
    yh = cell2mat(IndicatorVariablesInverse(YHatExact,Mod.Y.v));
    YHat{level}(which,:) = yh;
    
    [yn,~,yi] = unique(yh,'rows');
    cl = seekstring(ClasModel.ClassesName,yn);
    
    if level < Levels.NumLevels
        
        for i=1:size(yn,1)
            x = selectrow(X,yi==i);
            [YHat] = Classification(level+1,ClasModel.SubModel{cl(i)},...
                x,which(yi==i),YHat,Levels,PLSPar);
        end
    end
    
else
    if level < Levels.NumLevels
        YHat{level}(which,:) = repmat(ClasModel.ClassesName,length(which),1);
        [YHat] = Classification(level+1,ClasModel.SubModel{1},X,...
            which,YHat,Levels,PLSPar);
    else
        YHat{level}(which,:) = repmat(ClasModel.ClassesName,length(which),1);
    end
end