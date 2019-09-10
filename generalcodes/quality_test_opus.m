function [ZQuality]=quality_test_opus(ZSaisir)
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
%  First version: 12.03.05                                                             %
%  Revised: 24.01.10                                                                   %
%  Revised: 15.06.12                                                                   %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%  function [QTmatrix]=quality_test(ZSaisir);                                          %
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
%  Output:  
%                                                                                      %
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
ZSaisir1Deriv=saisir_derivative(ZSaisir,2,9,1);

        
% absorbance 1600-2100, here OPUS uses raw spectra
WN2=1600.0;
WN1=2100.0;
[y,i1]=min(abs(str2num(ZSaisir.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir.v)-WN2));
ZSaisir1600_2100=selectcol(ZSaisir,[i1:i2]);
[MaxAbs]=max(ZSaisir1600_2100.d,[],2);
[MinAbs]=min(ZSaisir1600_2100.d,[],2);

% (Signal 1) Amide 1: 1600-1700 (according to OPUS) - First derivative
WN2=1600.0;
WN1=1700.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir1600_1700=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxAmide1]=max(ZSaisir1600_1700.d,[],2);
[MinAmide1]=min(ZSaisir1600_1700.d,[],2);

% (Signal 2) Region 960-1200 - First derivative
WN2=960.0;
WN1=1200.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir960_1200=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxPoly]=max(ZSaisir960_1200.d,[],2);
[MinPoly]=min(ZSaisir960_1200.d,[],2);


% Water vapour: 1837-1847 (according to OPUS)- First derivative
WN2=1837.0;
WN1=1847.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir1837_1847=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxWater]=max(ZSaisir1837_1847.d,[],2);
[MinWater]=min(ZSaisir1837_1847.d,[],2);

% Noise: 2000-2100 (according to OPUS) - First derivative
WN2=2000.0;
WN1=2100.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir2000_2100=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxNoise]=max(ZSaisir2000_2100.d,[],2);
[MinNoise]=min(ZSaisir2000_2100.d,[],2);

% Ester: 1700-1800 (gives a measure for the fat content IN RAW SPECTRA)
WN2=1700.0;
WN1=1800.0;
[y,i1]=min(abs(str2num(ZSaisir.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir.v)-WN2));
ZSaisir1700_1800=selectcol(ZSaisir,[i1:i2]);
[MaxEster]=max(ZSaisir1700_1800.d,[],2);
[MinEster]=min(ZSaisir1700_1800.d,[],2);

% Fringes 2000-2200, here OPUS raw spectra
WN2=2000.0;
WN1=2300.0;
[y,i1]=min(abs(str2num(ZSaisir1Deriv.v)-WN1));
[y,i2]=min(abs(str2num(ZSaisir1Deriv.v)-WN2));
ZSaisir2000_2300=selectcol(ZSaisir1Deriv,[i1:i2]);
[MaxFringes]=max(ZSaisir2000_2300.d,[],2);
[MinFringes]=min(ZSaisir2000_2300.d,[],2);


Abs=MaxAbs-MinAbs;
Poly=MaxPoly-MinPoly;
Fringes=MaxFringes-MinFringes;
Noise=MaxNoise-MinNoise;
AmideI=MaxAmide1-MinAmide1;
AmideIN=AmideI./Noise;
PolyN=Poly./Noise;

IdentLowNoise= find(Noise<0.00000000001);
AmideIN(IdentLowNoise)=0.0;
PolyN(IdentLowNoise)=0.0;

Water=MaxWater-MinWater;
AmideW=AmideI./Water;
PolyW=(MaxPoly-MinPoly)./Water;

IdentLowWater=find(Water<0.00000000001);
AmideW(IdentLowWater)=0.0;
PolyW(IdentLowWater)=0.0;
Ester=MaxEster-MinEster;

[Nx Ny]=size(ZSaisir.d);
ZQuality.d=zeros(Nx,9);
ZQuality.d(:,1)=Abs;
ZQuality.d(:,2)=AmideI;
ZQuality.d(:,3)=Poly;
ZQuality.d(:,4)=Noise;
ZQuality.d(:,5)=AmideIN;
ZQuality.d(:,6)=PolyN;
ZQuality.d(:,7)=Water;
ZQuality.d(:,8)=AmideW;
ZQuality.d(:,9)=PolyW;
ZQuality.d(:,10)=Ester;
ZQuality.d(:,11)=Fringes;

ZQuality.i=ZSaisir.i;
ZQuality.v=['Absorbance '
            'AmideI     '
            'Poly       '
            'Noise      '
            'AmideIN    '
            'PolyN      '
            'Water Vapor'
            'AmideW     '
            'PolyW      '
            'Ester      '
            'Fringes    '];
          
          
          
          
          
          
          
          
          
          
          
