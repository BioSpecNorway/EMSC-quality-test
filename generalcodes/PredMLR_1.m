function Yhat = PredMLR_1(Xnew, mX, mY, beta)
% Computation of fitted values based on the q regresion models corresponding to the
% columns of <meanX>, <meanY> (mean values from traning data) and <beta> for new data Xnew:
[n, ~] = size(Xnew);
[~, k, q] = size(beta);

XnewS = Sentrer(Xnew,mX);
% Xnew   = Xnew - ones(n,1)*mX; % Centering of new data according to mean of training data
%
 k1 = size(mY,2);
Yhat = repmat(mY,n,k/k1,q);
Yhat = Yhat + mmat(XnewS,beta);
end