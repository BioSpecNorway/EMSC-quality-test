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
%  PLSR model for moulds data using OPUS quality test                                  %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Set path to ako's functions
Mainpath=genpath('C:\Users\valeta\Documents\codes-for-emsc-quality-test-paper-by-v.tafintseva\generalcodes');
addpath(Mainpath);
loadPath='\\nmbu.no\my\home\Documents\My Articles\EMSC QT paper\data\';
%%
DataInFile = 'MouldsNewName.mat';
load(strcat(loadPath,DataInFile))

ZX = delete_from_identifier(ZX,11,'cir');
ZX = delete_from_identifier(ZX,11,'her');
ZX = delete_from_identifier(ZX,11,'hem');
ZX = delete_from_identifier(ZX,11,'var');
ZX = delete_from_identifier(ZX,11,'plu');

ZX4Cal = delete_from_identifier(ZX,16,'04');
ZX4Cal = delete_from_identifier(ZX4Cal,16,'06');
ZX4Val1 = select_from_identifier(ZX,16,'04');
ZX4Val2 = select_from_identifier(ZX,16,'06');
ZX4Val = appendrow(ZX4Val1,ZX4Val2);

%% Opus quality test
[ZX4Cal_quality]=run_quality_test_opus(ZX4Cal);
[ZX4Val_quality]=run_quality_test_opus(ZX4Val);

%% ZXCal SPECTRA PRE-PROCESSING
% Pre-process spectra: derivative
ZXCal_2ndDer = saisir_derivative(ZX4Cal_quality,3,9,1); 

% SELECT RELEVANT REGIONS
[y,i1]=min(abs(str2num(ZXCal_2ndDer.v)-3050.0));
[y,i2]=min(abs(str2num(ZXCal_2ndDer.v)-2800.0));
ZX1=selectcol(ZXCal_2ndDer,[i1:i2]);
% 
[y,i1]=min(abs(str2num(ZXCal_2ndDer.v)-1800.0));
[y,i2]=min(abs(str2num(ZXCal_2ndDer.v)-900.0));
ZX2=selectcol(ZXCal_2ndDer,[i1:i2]);
ZXCalSelected=appendcol(ZX1,ZX2);

% Pre-process spectra: EMSC
[EMSCModel]=make_emsc_modfunc(ZXCalSelected);
[ZXCalEMSCCor,~,~]=cal_emsc(ZXCalSelected,EMSCModel,[],[]);
ZXCal = ZXCalEMSCCor;

%% ZXVal SPECTRA PRE-PROCESSING
% Pre-process spectra: derivative
ZXVal_2ndDer = saisir_derivative(ZX4Val_quality,3,9,1);

% SELECT RELEVANT REGIONS
[y,i1]=min(abs(str2num(ZXVal_2ndDer.v)-3050.0));
[y,i2]=min(abs(str2num(ZXVal_2ndDer.v)-2800.0));
ZX1=selectcol(ZXVal_2ndDer,[i1:i2]);
% 
[y,i1]=min(abs(str2num(ZXVal_2ndDer.v)-1800.0));
[y,i2]=min(abs(str2num(ZXVal_2ndDer.v)-900.0));
ZX2=selectcol(ZXVal_2ndDer,[i1:i2]);
ZXValSelected=appendcol(ZX1,ZX2);

[ZXValEMSCCor,~,~]=cal_emsc(ZXValSelected,EMSCModel,[],[]);
ZXVal = ZXValEMSCCor;

%% CLASSIFICATION

%% Tree specification
Levels.Names = ['Phylom ';'Class  ';'Family ';'Genus  ';'Species'];% 
Levels.NameRange{1,1} = 1:2;
Levels.NameRange{2,1} = 3:4;
Levels.NameRange{3,1} = 5:7;
Levels.NameRange{4,1} = 8:10;
Levels.NameRange{5,1} = 11:13;
Levels.NumLevels = size(Levels.Names,1);
if size(Levels.Names,1)~=size(Levels.NameRange,1)
    error('Levels information mismatch. Check levels names.')
end

N = size(ZXCal.d,1);
ZY = split2levels(Levels.NameRange,ZXCal);

%% Define parameters
% Define CV type:  'Full CV','Contiguous Blocks','Random Subsets',
% 'Venetian Blinds',
CVPar.CVType = 'Specified';%
CVPar.cvfold = 5;
CVPar.cvind = [16;17];
CVPar = CVParGenerator(CVPar,ZXCal);

%% Define PLS method
PLSPar.PLSMethod = 'pls';
PLSPar.PLSType = 'DA'; % Discriminant Analysis
PLSPar.p_critical = 0.05;     % Significance level for sensitivity in the model selection: measured in precent of optimal error
PLSPar.pc = 15;               % Maximum number of PC's

%% Training
level = 1;
which = (1:N)';
PlotTrain.ConfM = 0;
PlotTrain.RegrCoef = 0;
[ClasModel] = EstimateClassTree(level,ZXCal,ZY,which,PLSPar,CVPar,Levels,PlotTrain);

%% VALIDATION based on external set
Plot.Levels = 1;
Plot.Models = 0;
[ZYValHat,Results] = ClassEvaluation(ClasModel,ZXVal,Levels,PLSPar,Plot);
