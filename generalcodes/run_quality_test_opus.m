function [ZSaisirQT,IdentAll]=run_quality_test_opus(ZSaisir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      %
%                                                                                      %
%  Achim Kohler                                                                        %
%  Department of Mathematical Sciences and Technology (IMT)                            %
%  Norwegian University of Life Sciences                                               %
%  1432 Ås                                                                             %
%  Norway                                                                              %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%  First version: 20.06.16                                                             %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%  function [QTmatrix]=run_quality_test(ZSaisir,????);                                 %
%                                                                                      %
%  Runs a quality test on a Saisir structure                                           %
%  Literature: See Bruker user manual for microorganisms                               % 
%                                                                                      %
%                                                                                      %
%  Status: Running                                                                     %
%                                                                                      %
%  Input:   Saisir structure                                                           %
%           option (0: work with original spectra , 1: first derivative)               %                   
%                                                                                      %
%  Output:  Quality analysis                                                           %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate quality parameters
[ZQuality]=quality_test_opus(ZSaisir);

Abs=ZQuality.d(:,1);
AmideI=ZQuality.d(:,2);
Poly=ZQuality.d(:,3);
Noise=ZQuality.d(:,4);
AmideIN=ZQuality.d(:,5);
PolyN=ZQuality.d(:,6);
Water=ZQuality.d(:,7);
AmideIW=ZQuality.d(:,8);
PolyW=ZQuality.d(:,9);
% Ester=ZQuality.d(:,10);
Fringes=ZQuality.d(:,11);
% CHstretch=ZQuality.d(:,12);


N=size(ZSaisir.d,1);
IndicatorAll=zeros(N,1);


% Remove low absorbance spectra
IdentLowAbs= find(Abs<0.3450); %!!!!!!!!!!!! 0.3450
IndicatorAll(IdentLowAbs)=1;

% Remove high absorbance spectra
IdentHighAbs= find(Abs>1.2450); %!!!!!!!!!!!! 1.2450
IndicatorAll(IdentHighAbs)=1;

% Remove noisy spectra
IdentHighNoise= find(Noise>0.00015); %!!!!!!!!!!!!0.00015
IndicatorAll(IdentHighNoise)=1;

% Remove water vapor spectra
IdentHighWater= find(Water>0.0003);%!!!!!!!!!!!!0.0003
IndicatorAll(IdentHighWater)=1;

% Remove low AmideIN
IdentLowAmideIN= find(AmideIN<200);%!!!!!!!!!!!!200
IndicatorAll(IdentLowAmideIN)=1;

% Remove low PolyN
IdentLowPolyN= find(PolyN<40.0);%!!!!!!!!!!!!40
IndicatorAll(IdentLowPolyN)=1;

% Remove low AmideIW
IdentLowAmideIW= find(AmideIW<100.0);%!!!!!!!!!!!!100
IndicatorAll(IdentLowAmideIW)=1;

% Remove low PolyW
IdentLowPolyW= find(PolyW<20.0);%!!!!!!!!!!!!20
IndicatorAll(IdentLowPolyW)=1;

% % Remove high Fringes
% IdentHighFringes= find(Fringes>0.00003);%!!!!!!!!!!!!0.00003
% IndicatorAll(IdentHighFringes)=1;

IdentAll=find(IndicatorAll==1);

ZSaisirQT=deleterow(ZSaisir,IdentAll);






          
          
          
          
          
          
          
          
          
          
