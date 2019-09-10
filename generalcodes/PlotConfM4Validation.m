function PlotConfM4Validation(level,ClasModel,which,Levels,Y,YHat)
% function PlotConfM4Validation(level,ClasModel,which,Levels,Y,YHat)
% plot confusion matrices for each model when validating
% inputs:   level       number of the current level
%           ClasModel   hierarhical tree
%           which       sample index on the current level
%           Levels      structure containing information about the tree
%           Y           true values (char)
%           YHat        predictions (char)
%% V.Tafintseva 11.05.2016            
% 

NY = size(Y,1);
Nl = length(Levels.NameRange{level});

if ~isempty(ClasModel.Model)
    y = Y(:,Levels.NameRange{level});
    yhat = YHat(:,Levels.NameRange{level});
    
    Result = ClassResult(y,yhat);
  
    plotConfusionMat(Result{1})
    ax = axis;
    ty = text(ax(2),ax(3),['AOpt=',num2str(ClasModel.Model.AOpt)]);
    set(ty,'HorizontalAlignment','right','VerticalAlignment','top', ...
      'Rotation',90,'FontSize',12,'FontWeight','bold');
    set(gca,'FontSize',12,'FontWeight','bold')
    tx = text(ax(2),ax(4),'Val');
    set(tx,'HorizontalAlignment','right','VerticalAlignment','top', ...
      'FontSize',12,'FontWeight','bold');
 
  %%
  if level < Levels.NumLevels
  [yn,~,yi] = unique(y,'rows');
  nLev = size(yn,1);
  yj = strcmp(mat2cell(y,ones(NY,1),Nl),mat2cell(yhat,ones(NY,1),Nl));
   
    for i=1:nLev
     yf = find((yi==i) & (yj==1));
     if yf
       PlotConfM4Validation(level+1,ClasModel.SubModel{i},which(yf),...
                          Levels,Y(yf,:),YHat(yf,:));
     end
    end
  end
   
else 
  if level < Levels.NumLevels
   % YHat(which,:) = repmat(ClasModel.ClassesName,length(which),1);
    PlotConfM4ValidationDA(level+1,ClasModel.SubModel{1},which,Levels,Y,YHat);
  end

end 

end