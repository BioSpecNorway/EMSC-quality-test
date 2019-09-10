function [ZImageCorrected,ZImageResiduals,ZImageEMSCParam]=cal_emscImage(ZImage,EMSCModel,Mask,option,ZWeights);
%cal_emscImage  [ZImageCorrected,ZImageResiduals,ZImageEMSCParam]=cal_emscImage(ZImage,EMSCModel,Mask,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      %
%  Achim Kohler                                                                        %
%                                                                                      %
%  Biospectroscopy and Data Modelling Group                                            %
%  Dept. of Mathematical Sciences and Technology                                       %
%  Norwegian University of Life Sciences (www.umb.no)                                  %
%                                                                                      %
%  Homepage: http://arken.umb.no/~achik/                                               %
%                                                                                      %
%  First version: 15.06.2012                                                           %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%  function [ZOutput]=example_func(ZInput);                                            %
%                                                                                      %
%  option=0: only correction of additive terms                                         %
%  If option is not defined or 1: default is EMSC with correction of additive and      %
%                                 multiplicative variations                            %   
%                                                                                      %
%  Literature: See ......                                                              %
%                                                                                      %
%                                                                                      %
%  Status: Running                                                                     %
%                                                                                      %
%  Input:                                                                              %
%                                                                                      %
%                                                                                      %
%  Output:                                                                             %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[N,M,LSpec]=size(ZImage.d);

% check the input
if (isempty(option))
    %default: all model function
    Option=1;
else
    Option=option;
end


if (isempty(Mask))
    [N,M,L]=size(ZImage.d);
    ZSaisir.d=reshape(ZImage.d,N*M,L);
    ZSaisir.v=ZImage.v;
    ZSaisir.i=num2str((1:N*M)');

    [ZCorrected,ZResiduals,ZEMSCParam]=cal_emsc(ZSaisir,EMSCModel,ZWeights,option);

    ZImageCorrected=ZImage;
    ZImageCorrected.d=reshape(ZCorrected.d,N,M,L);
    ZImageCorrected.v=ZCorrected.v;

    ZImageResiduals=ZImage;
    ZImageResiduals.d=reshape(ZResiduals.d,N,M,L);
    ZImageResiduals.v=ZCorrected.v;

    [L]=size(ZEMSCParam.d,2);
    ZImageEMSCParam=ZImage;
    ZImageEMSCParam.d=reshape(ZEMSCParam.d,N,M,L);
    ZImageEMSCParam.v=ZEMSCParam.v;


else
    
    
    [Nx,Ny,L]=size(ZImage.d);
    N=Nx*Ny;
    MaskVector=reshape(Mask,[N,1]);
    ZSaisirVector.d=reshape(ZImage.d,N,L);
    ZSaisirVector.i=num2str((1:N)');
    ZSaisirVector.v=ZImage.v;
    
    [Index]=find(MaskVector);
    ZSaisirRed.d(1:size(Index,1),:)=ZSaisirVector.d(Index,:);
    ZSaisirRed.v=ZImage.v;
    ZSaisirRed.i=num2str((1:size(Index,1))');
        
   [ZCorrected,ZResiduals,ZEMSCParam]=cal_emsc(ZSaisirRed,EMSCModel,ZWeights,option);
   
   % ZImageCorrected
   ZSaisirVector.d=zeros(N,L);
   ZSaisirVector.d(Index,:)=ZCorrected.d(1:size(Index,1),:);
   ZImageCorrected.x=ZImage.x;
   ZImageCorrected.y=ZImage.y;
   ZImageCorrected.v=ZCorrected.v;
   ZImageCorrected.d=reshape(ZSaisirVector.d,Nx,Ny,L);
   
   % ZImageResiduals
   ZSaisirVector.d=zeros(N,L);
   ZSaisirVector.d(Index,:)=ZResiduals.d(1:size(Index,1),:);
   ZImageResiduals.x=ZImage.x;
   ZImageResiduals.y=ZImage.y;
   ZImageResiduals.v=ZCorrected.v;
   ZImageResiduals.d=reshape(ZSaisirVector.d,Nx,Ny,L);
   
   % ZImageParam
   [LPar]=size(ZEMSCParam.d,2);
   ZSaisirVector.d=zeros(N,LPar);
   ZSaisirVector.d(Index,:)=ZEMSCParam.d(1:size(Index,1),:);
   ZImageEMSCParam.x=ZImage.x;
   ZImageEMSCParam.y=ZImage.y;
   ZImageEMSCParam.v=ZEMSCParam.v;
   ZImageEMSCParam.d=reshape(ZSaisirVector.d,Nx,Ny,LPar);
   
 
end



