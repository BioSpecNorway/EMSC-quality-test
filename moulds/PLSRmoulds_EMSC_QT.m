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
%  PLSR model for moulds data using EMSC quality test                                  %
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

%% EMSC Quality control
minb = 0.345;
maxb = 1.245;

[~,index_ZXqt] = EMSC_quality_test(ZX4Cal,minb,maxb);
ZX4Cal_quality = deleterow(ZX4Cal,index_ZXqt);

%% ZXCal SPECTRA PRE-PROCESSING
% Pre-process spectra: derivative
ZXCal_2ndDer = saisir_derivative(ZX4Cal_quality,3,9,1);

% SELECT RELEVANT REGIONS
[y,i1]=min(abs(str2num(ZXCal_2ndDer.v)-3050.0));
[y,i2]=min(abs(str2num(ZXCal_2ndDer.v)-2800.0));
ZXCal1=selectcol(ZXCal_2ndDer,[i1:i2]);
%
[y,i1]=min(abs(str2num(ZXCal_2ndDer.v)-1800.0));
[y,i2]=min(abs(str2num(ZXCal_2ndDer.v)-900.0));
ZXCal2=selectcol(ZXCal_2ndDer,[i1:i2]);

ZXCalSelected=appendcol(ZXCal1,ZXCal2);

[EMSCModel]=make_emsc_modfunc(ZXCalSelected);
[ZXCalEMSCCor,~,~]=cal_emsc(ZXCalSelected,EMSCModel,[],[]);
ZXCal = ZXCalEMSCCor;

%% ZXVal SPECTRA PRE-PROCESSING
[ZX4Val_quality]=run_quality_test_opus(ZX4Val);
% Pre-process spectra: derivative
ZXVal_2ndDer = saisir_derivative(ZX4Val_quality,3,9,1);

% SELECT RELEVANT REGIONS
[y,i1]=min(abs(str2num(ZXVal_2ndDer.v)-3050.0));
[y,i2]=min(abs(str2num(ZXVal_2ndDer.v)-2800.0));
ZXVal1=selectcol(ZXVal_2ndDer,[i1:i2]);
%
[y,i1]=min(abs(str2num(ZXVal_2ndDer.v)-1800.0));
[y,i2]=min(abs(str2num(ZXVal_2ndDer.v)-900.0));
ZXVal2=selectcol(ZXVal_2ndDer,[i1:i2]);
ZXValSelected=appendcol(ZXVal1,ZXVal2);

% Pre-process spectra: EMSC
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
CVPar.cvind = [16,17];
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

%%
Plot.Levels = 1;
Plot.Models = 0;
[ZYValHat,Results] = ClassEvaluation(ClasModel,ZXVal,Levels,PLSPar,Plot);

function [ZX_EMSC_quality,IdentAll] = EMSC_quality_test(ZX,minb,maxb)

N = size(ZX.d,1);
IndicatorAll = zeros(N,1);

Reference = correctReference(ZX);
[EMSCModel] = make_emsc_modfunc(Reference);
[ZXEMSCCor,~,EMSCParam] = cal_emsc(ZX,EMSCModel,[],[]);

IdentLow_b = find(EMSCParam.d(:,4)<minb);%
IndicatorAll(IdentLow_b)=1;

IdentHigh_b = find(EMSCParam.d(:,4)>maxb);%
IndicatorAll(IdentHigh_b)=1;

IdentAll = find(IndicatorAll==1);
ZX_EMSC_quality = deleterow(ZXEMSCCor,IdentAll);
end

