function Y = IndicatorVariablesInverse(YhatExact,groupNames)
% function Y = IndicatorVariablesInverse(YhatExact,groupNames)
% tranforms exact Yhat to string group label
% inputs:   YhatExact       exactly estimated prediction values 
%           groupNames      vector of names of groups
% outputs:  Y               1xpc cell of string-labelled predictions  
%                           for each principal component
%% V.Tafintseva 11.05.2016
%
[NY,KY,pc] = size(YhatExact);
Y = cell(1,pc);

%% Rearranging YhExact
IndicatorMatTemp = num2cell(YhatExact,[1 2]); 
YTemp = squeeze(IndicatorMatTemp); 

for j = 1:pc
  [~,index] = min(abs(ones(NY,KY)-YTemp{j}),[],2); 
  Y{j} = groupNames(index,:);
end

end

