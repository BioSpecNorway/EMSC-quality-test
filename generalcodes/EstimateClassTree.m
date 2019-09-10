function [Model] = EstimateClassTree(level,ZX,Yprim,which,PLSPar,CVPar,Levels,Plot)
% function [Model] = EstimateClassTree(level,ZX,Yprim,Yadd,which,PLSPar,CVPar,Levels,Plot)
% recursive function builds the hierarchical tree according to Levels
% inputs:   level       current level number
%           ZX          saisir data
%           Yprim       1xNumLevels cell with sample names       
%           which       sample indeces used for calibration
%           PLSPar      parameters of PLS regression
%           CVPar       parameters of Cross Validation
%           Levels      hierarhical tree info (Names, NameRange, NumLevels)
%           Plot        structure (ConfM and RegrCoef)
% outputs:  Model       hierarhical classification tree
%   
%% V.Tafintseva 11.05.2016

orig = unique(Yprim{level},'rows');
if size(orig,1) > 1
    y.i =  Yprim{level};
    ZY = IndicatorVariables(y);
    PLSPar.ClasLevel = Levels.Names(level,:);
    ClassModel = CV_PLSDA(ZX,ZY,PLSPar.pc,PLSPar,CVPar,Plot);
else
end

if level < Levels.NumLevels % Split into submodels recursively
    [C,~,ic] = unique(Yprim{level},'rows');
    nLev = size(C,1);
    if nLev > 1
        Model.ClassesName = orig;
        Model.Model = ClassModel;
        Model.SubModel = cell(nLev,1);
    else
        Model.ClassesName = Yprim{level}(1,:);
        Model.Model = [];
        Model.SubModel = [];
    end
    for i=1:nLev
        x = selectrow(ZX,ic==i);
        yprim = Yprim;
        for j = 1:length(Yprim)
            yprim{j} = yprim{j}(ic==i,:);
        end
        
        cvPar = CVPar;
        cvPar.cv = CVPar.cv(ic==i); %[~,~,cv] = unique(Cv);
        if size(unique(yprim{level+1},'rows'),1) > 1 || level < length(Yprim)
            [SubModel] = EstimateClassTree(level+1,x,yprim,...
                which(ic==i),PLSPar,cvPar,Levels,Plot);
            Model.SubModel{i} = SubModel;
        else
            Model.SubModel{i} = [];
        end
    end
else
    if size(unique(Yprim{level},'rows'),1) > 1
        Model.ClassesName = orig;
        Model.Model = ClassModel;
    else
        Model.ClassesName = Yprim{level}(1,:);
        Model.Model = [];
    end
end

end