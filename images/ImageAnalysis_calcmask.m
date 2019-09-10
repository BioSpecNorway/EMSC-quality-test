%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      %
%  Valeria Tafintseva                                                                  %
%  Norwegian University of Life Sciences                                               %
%                                                                                      %
%  codes for article "Extended  multiplicative signal correction for quality control   %
%                       of FTIR spectra and images of biological material"             %
%  09.09.19                                                                            %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%  Image analysis, finding segmentation (masks) using EMSC models                      %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Set path to ako's functions
Mainpath=genpath('\\nmbu.no\my\home\Documents\My Articles\EMSC QT paper\codes2publish\generalcodes');
addpath(Mainpath);
loadPath='\\nmbu.no\my\home\Documents\My Articles\EMSC QT paper\data\images\';
%%
Filename='FPA_17_01_05_Muccir_phosph1_day2_54h_repl2_RAWD.ir5';
% Return information about HDF5 file
info=h5info(strcat(loadPath,Filename));

%% Build a saisir structure for each projection
ZProjection.d = squeeze(h5read(strcat(loadPath,Filename),...
    '/datacube',[1,1,1,1],[128,128,1,765]));

% Read the wavenumbers
ZProjection.v = num2str(h5read(strcat(loadPath,Filename),'/wavenumbers'));
[Nx,Ny,LSpec] = size(ZProjection.d);
ZProjection.x = num2str((1:Nx)');
ZProjection.y = num2str((1:Ny)');

[~,i2] = min(abs(str2num(ZProjection.v)-3010));
[~,i1] = min(abs(str2num(ZProjection.v)-3010));
Zsection.d = ZProjection.d(:,:,i1:i2);
Zsection.i = '3010';

figure;
pcolor(Zsection.d);
colorbar;
set(gca,'FontSize',12,'FontWeight','bold','linewidth',1)
title([Zsection.i, 'cm^{-1}']);

% Build the EMSC model and correct baseline
ModelOption = 1;
[EMSCModel] = make_emsc_modfuncImage(ZProjection,ModelOption,[]);

%% Calculate the EMSC model parameters
% Define Weights: 1 for relevant regions, 0.0001 for unimportant ones
[~, M, L] = size(ZProjection.d);
ZWeights.d = ones(1,L)*0.0001;
ZWeights.i = 'Weights';
ZWeights.v = ZProjection.v;

[~,i2] = min(abs(str2num(ZWeights.v)-898.64));
[~,i1] = min(abs(str2num(ZWeights.v)-2000.0));
ZWeights.d(1,i1:i2) = 1.0;

[~,i2] = min(abs(str2num(ZWeights.v)-2600.0)); 
[~,i1] = min(abs(str2num(ZWeights.v)-3845.28)); 
ZWeights.d(1,i1:i2) = 1.0;
    
[~,ZImageResiduals,ZImageEMSCParam] = cal_emscImage(ZProjection,...
    EMSCModel,[],[],ZWeights);

%% Plot the residual image
[Nx,Ny,Nz] = size(ZImageResiduals.d);
N = Nx*Ny;
ImageVector = reshape(ZImageResiduals.d,[N,Nz]);
EMSC_RMSE_vector = sqrt(sum(ImageVector.^2,2)/Nz);
EMSC_RMSE_image = squeeze(reshape(EMSC_RMSE_vector,[Nx,Ny,1]));

figure;
pcolor(EMSC_RMSE_image);
colorbar;
set(gca,'FontSize',12,'FontWeight','bold','linewidth',1)
title('RMSE, residual after first EMSC');
    
%% Plot the EMSC model parameter images
figure;
subplot(2,2,1)
ImageA = squeeze(ZImageEMSCParam.d(:,:,1));
pcolor(ImageA);
colorbar;
set(gca,'FontSize',12,'FontWeight','bold','linewidth',1)
title('baseline');
subplot(2,2,2)
ImageD = squeeze(ZImageEMSCParam.d(:,:,2));
pcolor(ImageD);
colorbar;
set(gca,'FontSize',12,'FontWeight','bold','linewidth',1)
title('linear');
subplot(2,2,3)
ImageE = squeeze(ZImageEMSCParam.d(:,:,3));
pcolor(ImageE);
colorbar;
set(gca,'FontSize',12,'FontWeight','bold','linewidth',1)
title('quadratic');
subplot(2,2,4)
ImageB = squeeze(ZImageEMSCParam.d(:,:,end)); %
pcolor(ImageB);
colorbar;
set(gca,'FontSize',12,'FontWeight','bold','linewidth',1)
title('multiplicative');
suptitle('EMSC model parameters')

%% Develop mask on multiplicative effect
ImageB = squeeze(ZImageEMSCParam.d(:,:,4));
[Nx,Ny] = size(ImageB);
N=Nx*Ny;
ImageVector = reshape(ImageB,[N,1]);

b=2.2;

EMSCMask = calcmask(ImageB,b);
figure
imshow(EMSCMask,'InitialMagnification', 400);
set(gca,'FontSize',10,'FontWeight','bold','linewidth',1,'YDir','normal');
imwrite(EMSCMask,['maskb',num2str(b),'.png'])
title(['Mask using b=',num2str(b)])

%%
function mask = calcmask(ImageB,par)
 
 [Nx,Ny] = size(ImageB);
 N=Nx*Ny;
 ImageVector = reshape(ImageB,[N,1]);
 ImageMaskVector = ones(N,1);
 Indmask = find(ImageVector<par);   %% to be set % Threshold needs to be optimized
 ImageMaskVector(Indmask) = 0.0; % Set the mask
 mask = reshape(ImageMaskVector,[Nx,Ny]); % Fold the mask
 
 end