function [yModel] = CV_PLSDA(X,Yprim,k,PLSPar,CVPar,Plot)
% function [yModel] = CV_MBPLSDA(X,Yprim,Yadd,k,PLSPar,CVPar, Plot)
% Purpose: Establish a calibration model using CV for Y indicator matrix
% 
% Input:   X       saisir structure
%          Yprim   cell containing names of samples on the PLSPar.ClasLevel
%          k       number of PCs
%          PLSPar  parameters of PLS regression
%          CVPar   parameters of Cross Validation
%          Plot    switch with Plot.RegrCoef and Plot.CoonfMat
%
% Output:  yModel      Calibration model for y
%  
%% V.Tafintseva 27.04.2015

[NX, KX] = size(X.d);
[NY, KY] = size(Yprim.d);
if NX ~= NY
    error('namber of samples in X and Y matrices should be same')
end
cvfold = CVPar.cvfold;
cv = CVPar.cv;
ClasLevel = PLSPar.ClasLevel;
p_critical = PLSPar.p_critical; % for error optimization % of opt error

% ------------------------------------------------------------------------
% ---------------- Loop for crossvalidation & full modeling ---------------
% -------------------------------------------------------------------------
YhExactcv = ones(NX,KY,k)*NaN; % Responses to be estimated by cross-validation
Ycv.i  = repmat(' ',NX,size(Yprim.i,2)); % Responses given
Ycv.d  = ones(NY,KY)*NaN;
Ycv.v  = Yprim.v;
MCRcv = zeros(1,k);

%% Cross-validation
for i = 1:cvfold  % cv-loop , cvfold - number of segments
    
    cvout = find(cv==i);                      % Objects taken out
    cvin  = setdiff((1:NX)',cvout);               % Objects held in
    % Building models from observations "held in":
    Xin = selectrow(X,cvin);
    Yprimin = selectrow(Yprim,cvin);
    [Beta,~,~,~,~,mXin, mYin] = pls2(Xin,Yprimin,k);
    
    YhExactcv(cvout,:,:) = PredMLR_1(X.d(cvout,:), mXin.d, mYin.d, Beta.d); % Predictions: precise estimation
    Ycv.i(cvout,:) = Yprim.i(cvout,:);  % and copy of corresponding input responses:
    Ycv.d(cvout,:) = Yprim.d(cvout,:);
    Ycv.f(cvout,:) = X.i(cvout,:);
end           % end cv-loop

%% Classification by CV: Cmat - confusion matrix, Yh - predictions 'char' vector
Yhcv = IndicatorVariablesInverse(YhExactcv,Yprim.v);

[ClasResultcv] = ClassResult(Ycv.i,Yhcv); %

% ERROR calculation
for j=1:k
    MCRcv(j) = 1-ClasResultcv{j}.statistics.accuracy;
end

[~,AMin] = find(MCRcv == min(MCRcv),1,'first');
cvsignificant = binocdf(NX*min(MCRcv),NX*ones(1,k), MCRcv);
AOpt = find(cvsignificant > p_critical,1,'first');

if Plot.ConfM == 1
    plotConfusionMat(ClasResultcv{AOpt},'cal')
    ax = axis;
    ty = text(ax(2),ax(3),['AOpt=',num2str(AOpt)]);
    set(ty,'HorizontalAlignment','right','VerticalAlignment','top', ...
        'Rotation',90,'FontSize',12,'FontWeight','bold');
    set(gca,'FontSize',12,'FontWeight','bold')
    tx = text(ax(2),ax(4),'Cal');
    set(tx,'HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',12,'FontWeight','bold');
    title(ClasLevel)
end

%% Final full modelling based on all observations:
[Beta, W, P, Q, T, mX, mY] = pls2(X, Yprim, k);

% ------- FINAL REGRESSION COEFFS ---------------------
if AOpt > 0
    B.d = Beta.d(:,:,AOpt);  % regression coefficients at AOpt
    B.i = X.v;
    B.v = Yprim.v;
    B0.d = mY.d - mX.d*B.d;
    B0.v = Yprim.v;
    B0.i = 'slope';
    MCRAOpt = MCRcv(AOpt);
    ClasResultAOpt = ClasResultcv{AOpt};
else
    B = zeros(KX,KY);
    B0 = mY;
    MCRAOpt = 0;
    ClasResultAOpt = [];
end

if AMin>0
    MCRAMin = MCRcv(AMin);
else
    MCRAMin = 0;
end % if AMin

if Plot.RegrCoef == 1
    figure
    set(gcf,'Color',[1 1 1]);
    plot(str2num(B.i),B.d,'LineWidth',1.5)
    set(gca,'XDir','reverse','FontSize',12,'LineWidth',1.5);
    axis tight;
    xlabel('Wavenumber [cm^-^1]','FontSize',14,'FontWeight','Bold');
    title(['Regression coefficients for ',ClasLevel],'FontSize',14);
    leggend = [];
    for j = 1:KY
        leggend = [leggend,'''',Yprim.v(j,:),'''',','];
    end
    leggend(end)=[];
    eval(['legend(',leggend,',''Location'',''Best'')'])
end



% Saving results:
yModel.ClasLevel=ClasLevel;
yModel.PLSPar=PLSPar;
yModel.Y=Yprim;
yModel.mX=mX;
yModel.mY=mY;
yModel.AMin=AMin;
yModel.AOpt=AOpt;
yModel.W=W;
yModel.P=P;
yModel.Q=Q;
yModel.T=T;
yModel.B=B;
yModel.B0=B0;
yModel.Yhcv=Yhcv(:,AOpt);
yModel.Ycv=Ycv;
yModel.ClasResultAOpt=ClasResultAOpt;
yModel.MCRcv=MCRcv;
yModel.MCRAOpt=MCRAOpt;
yModel.MCRAMin=MCRAMin;
end




