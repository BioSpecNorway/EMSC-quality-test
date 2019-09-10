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
%  RF model for yeasts data using OPUS quality test                                    %
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
DataInFile = 'SABfullName4levels';
load(strcat(loadPath,DataInFile))

ZX4Cal = delete_from_identifier(ZX,19,'01');
ZX4Cal = delete_from_identifier(ZX4Cal,19,'05');
ZX4Val1 = select_from_identifier(ZX,19,'01');
ZX4Val2 = select_from_identifier(ZX,19,'05');
ZX4Val = appendrow(ZX4Val1,ZX4Val2);

%% Opus quality test
[ZX4Cal_quality]=run_quality_test_opus(ZX4Cal);
[ZX4Val_quality]=run_quality_test_opus(ZX4Val);

%% ZXCal SPECTRA PRE-PROCESSING
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
% ZXSelected=ZX1;

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
ZXVal1=selectcol(ZXVal_2ndDer,[i1:i2]);
%
[y,i1]=min(abs(str2num(ZXVal_2ndDer.v)-1800.0));
[y,i2]=min(abs(str2num(ZXVal_2ndDer.v)-900.0));
ZXVal2=selectcol(ZXVal_2ndDer,[i1:i2]);

ZXValSelected=appendcol(ZXVal1,ZXVal2);
[ZXValEMSCCor,~,~]=cal_emsc(ZXValSelected,EMSCModel,[],[]);
ZXVal = ZXValEMSCCor;

%% Random forest
ind_groups = 9:11; N = size(ZXCal.d,1);
samplesnames = ZXCal.i(:,ind_groups);
names2RF = mat2cell(samplesnames,ones(1,N),size(samplesnames,2));

rng default
numtrees = 300;
b = TreeBagger(numtrees,ZXCal.d,names2RF,'oobpred','on');
OOBer = oobError(b);
OOBmargin = oobMargin(b);

%% Predictions
% Calibration result
YCal = ZXCal.i(:,ind_groups);
YCalhat = cell2mat(predict(b,ZXCal.d));

ClasResultCal = ClassResult(YCal,YCalhat);
plotConfusionMat(ClasResultCal,'cal')

% Validation result
YVal = ZXVal.i(:,ind_groups);
YValhat = cell2mat(predict(b,ZXVal.d));

ClasResultVal = ClassResult(YVal,YValhat);
plotConfusionMat(ClasResultVal,'val')