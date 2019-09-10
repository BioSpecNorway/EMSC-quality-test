function [res] = ClassResult(Y,Yhat)

if ~iscell(Yhat)
    Yh{1} = Yhat;
else Yh = Yhat;
end

[~,k] = size(Yh);

res = cell(1,k);

for j = 1:k
    if isempty(Yh{j})
        res{j}.Y = Y;
        res{j}.Yh = Yh{j};
        res{j}.statistics.ModelMCR = 10000;
        continue
    end
  statistics = confusionmatStats(Y,Yh{j});
  res{j}.Y = Y;
  res{j}.Yh = Yh{j};
  res{j}.statistics = statistics;
end
end