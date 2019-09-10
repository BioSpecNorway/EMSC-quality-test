function [beta, W, P, Q, T, mX, mY] = pls2(ZX,ZYprim,pc)
%basic_pls2            - PLS2 on several variables, several dimensions
%% V.Tafintseva 23.06.2016

X = ZX.d;
Yprim = ZYprim.d;

[NX, KX] = size(X);
[NY, KY] = size(Yprim);
% Yprim = Yprim./(ones(NY,1)*sqrt(sum(Yprim)));
mY_t = mean(Yprim);
mX_t = mean(X);

X = X-ones(NX,1)*mX_t;
Yprim=Yprim-ones(NX,1)*mY_t;  
% Declaration of variables
W_t = zeros(KX,pc);
T_t = zeros(NX,pc);
P_t = zeros(KX,pc);
Q_t = zeros(KY,pc);
beta_t=zeros(KX,KY,pc);

for i=1:pc
    XY=X'*Yprim;      
    [u,~,~]=svd(XY,0);    
    w=u(:,1);
    W_t(:,i)=w; % store the orthogonal loading weights
    
    t=X*w;
   	T_t(:,i)=t;	 % common scores  
    nt=t/(t'*t);  % norm of t
    q=Yprim'*nt;  % calculate the loadings of y
    Q_t(:,i)=q;  % store the loadings of y
    p=X'*nt;  % calculate the loadings of X
    P_t(:,i)=p; % store the loadings of X
    
    beta_t(:,:,i) = W_t(: ,1:i) * pinv(P_t(:,1:i)'*W_t(: ,1:i)) * Q_t(:,1:i)';
    
    X=X-t*p';  % deflation of X
    Yprim=Yprim-t*q';  % deflation of Y
    
end

mX.d = mX_t;
mX.v = ZX.v;
mX.i = 'mean X';

mY.d = mY_t;
mY.v = ZYprim.v;
mY.i = 'mean Y';

T.d = T_t;
T.i = ZX.i;
T.v = [char(ones(pc,1))*'PC ', num2str((1:pc)')];

W.d = W_t;
W.i = ZX.v;
W.v = [char(ones(pc,1))*'PC ', num2str((1:pc)')];

P.d = P_t;
P.i = ZX.v;
P.v = [char(ones(pc,1))*'PC ', num2str((1:pc)')];

Q.d = Q_t;
Q.i = ZYprim.v;
Q.v = [char(ones(pc,1))*'PC ', num2str((1:pc)')];

beta.d = beta_t;
beta.i = ZX.v;
beta.v = ZYprim.v;
beta.l = [char(ones(pc,1))*'PC ', num2str((1:pc)')];
end

