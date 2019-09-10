function [EMSC]=make_emsc_modfunc(Zsaisir,option);
%make_emsc_modfunc   [EMSC]=make_emsc_modfunc(Zsaisir,option)
%
%  Establishes the basic EMSC model without any extensions in addition to
%  polynomials
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      %
%  Achim Kohler                                                                        %
%  Center for Biospectroscopy and Data Modelling                                       %
%  Nofima Mat                                                                          %
%  Norwegian Food Research Institute                                                   %
%  Osloveien 1                                                                         %
%  1430 Ås                                                                             %
%  Norway                                                                              %
%                                                                                      %
%  12.03.05                                                                            %
%  18.03.09 short revision                                                             %                                                                                      %
%                                                                                      %
%  todo: normalise polynomial functions (not really necessary)                         %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%  function [EMSC]=make_emsc_modfunc(Zsaisir,option);                                  %
%                                                                                      %
%  option=0: only baseline (no MSC/EMSC)                                               %
%  If option is not defined or 1: default is EMSC with linear and quadratic effect     %
%  option=2: MSC plus linear                                                           %
%  option=3: MSC                                                                       %
%  option=4: EMSC + cubic                                                              % 
%  option=5: EMSC + cubic + fourth order                                               % 
%  option=6: EMSC + cubic + fourth order + fifth order (not defined)                   % 
%  option=7: EMSC + cubic + fourth order + fifth order + sixth order                   %
%                                                                                      %
%  MeanScaling: Scale the Mean spectrum to the same size as the other basic EMSC par.  %                                                                                     %
%                                                                                      %
%                                                                                      %
%  Creates the basic emsc modell fucntions:                                            %
%  baseline,linear,quadratic,reference (in this oder)                                  %
%                                                                                      %
%                                                                                      %
%  Input:   Zsaisir                                                                    %
%                                                                                      %
%  Output:  EMSC structure with modell functions (see def. at end of file)             % 
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanScaling=0;


% check the input
if (nargin==1)
    %default: all model function
    Option=1;
elseif (nargin==2)
    Option=option;
end


%% Calculate the basic model functions
[Nx Ny]=size(Zsaisir.d);
WaveNum=str2num(Zsaisir.v);
Start=WaveNum(1);
End=WaveNum(Ny);

C=0.5*(Start+End);
M0=2.0/(Start-End);
M=4.0/((Start-End)*(Start-End));

for i=1:Ny
   Baseline(i)=1.0;
   Linear(i)=M0*(Start-WaveNum(i))-1.0;
   Quadratic(i)=M*(WaveNum(i)-C)*(WaveNum(i)-C);
   Cubic(i)=M*(1/(Start-End))*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C); %% check norm
   FourthOrder(i)=M*M*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C); %% check norm
   FifthOrder(i)=M*M*M0*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C); %% check norm
   SixthOrder(i)=M*M*M*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C)*(WaveNum(i)-C); %% check norm
   
end
Mean=mean(Zsaisir.d,1); % mean spectra

%% If necessary the mean spectrum can be scaled here! 
if (MeanScaling)
    MaxVal=max(Mean);
    MinVal=min(Mean);
    Mean=Mean/(MaxVal-MinVal);
end

%% collect the basic model functions
if (Option==0)
    MModel=[Baseline];
    ModelSpecNames=['Baseline        '];
elseif (Option==1)
    MModel=[Baseline;Linear;Quadratic;Mean];
    ModelSpecNames=['Baseline        ',
                    'Linear          ',
                    'Quadratic       ',
                    'Reference       '];
elseif (Option==2)
    MModel=[Baseline;Linear;Mean];
    ModelSpecNames=['Baseline        ',
                    'Linear          ',
                    'Reference       '];
elseif (Option==3)
    MModel=[Baseline;Mean]; % MSC
    ModelSpecNames=['Baseline        ',
                    'Reference       '];
elseif (Option==4)
    MModel=[Baseline;Linear;Quadratic;Cubic;Mean];
    ModelSpecNames=['Baseline        ',
                    'Linear          ',
                    'Quadratic       ',
                    'Cubic           ',
                    'Reference       '];
elseif (Option==5)
    MModel=[Baseline;Linear;Quadratic;Cubic;FourthOrder;Mean];
    ModelSpecNames=['Baseline        ',
                    'Linear          ',
                    'Quadratic       ',
                    'Cubic           ',
                    'fourth order    ',
                    'Reference       '];
elseif (Option==7)
    MModel=[Baseline;Linear;Quadratic;Cubic;FourthOrder;FifthOrder;SixthOrder;Mean];
    ModelSpecNames=['Baseline        ',
                    'Linear          ',
                    'Quadratic       ',
                    'Cubic           ',
                    'fourth order    ',
                    'fifth order     ',                    
                    'sixth order     ',                                        
                    'Reference       '];                
               
end
[L,O]=size(MModel);

EMSC.Model=MModel';
EMSC.ModelSpecNames=ModelSpecNames;
EMSC.NumBasicModelFunc=L;  % this defines the order of basic mod func
EMSC.NgoodSpec=0;          % Number of good spectra  
EMSC.NbadSpec=0;           % Number of bad spectra in Modell function
EMSC.NModelFunc=L;         % total number of model functions
EMSC.ModelVariables=Zsaisir.v;

