function [EMSCModel]=make_emsc_modfuncImage(ZsaisirImage,option,Mask);
%make_emsc_modfuncImage   [EMSC]=make_emsc_modfuncImage(ZsaisirImage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      %
%  Achim Kohler                                                                        %
%  Department of Mathematical Sciences and Technology (IMT)                            %
%  Norwegian University of Life Sciences                                               %
%  1432 Ås                                                                             %
%  Norway                                                                              %
%                                                                                      %
%                                                                                      %
%  26.10.2016                                                                          %
%  revised: 27.10.2016 extended to use option                                          %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%                                                                                      %
%  Description                                                                         %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%  Input:                                                                              %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%  Output:                                                                             %
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check the input
if isempty(option)
    %default: all model function
    Option=1;
else
    Option=option;
end

if isempty(Mask)
  
    
    [Nx,Ny,L]=size(ZsaisirImage.d);
    N=Nx*Ny;
    ZSaisir.d=reshape(ZsaisirImage.d,N,L);
    ZSaisir.v=ZsaisirImage.v;
    ZSaisir.i=num2str((1:N)');
    
    [EMSCModel]=make_emsc_modfunc(ZSaisir,Option);
    
    
    
else
    
    [Nx,Ny,L]=size(ZsaisirImage.d);
    N=Nx*Ny;
    MaskVector=reshape(Mask,[N,1]);
    ZSaisir.d=reshape(ZsaisirImage.d,N,L);
    
    [Index]=find(MaskVector);
    ZSaisirRed.d(1:size(Index,1),:)=ZSaisir.d(Index,:);
    ZSaisirRed.v=ZsaisirImage.v;
    ZSaisirRed.i=num2str((1:size(Index,1))');
    
    [EMSCModel]=make_emsc_modfunc(ZSaisirRed,Option);
    
end






