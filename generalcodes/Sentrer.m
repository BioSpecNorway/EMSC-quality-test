function [X, n, p] = Sentrer(X,mX)
[n,p] = size(X);
% Sentrering av X, samme som: %X = X-ones(n,1)*mX;
for i = 1:p
    X(:,i)=X(:,i)-mX(i);
end