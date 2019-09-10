function [refcorr,peak] = correctReference(saisir)
% function saisircorr = correct4b(saisir)
% made to correct for differences in spectra prior to quality test
% with EMSC parameter b
[~,K] = size(saisir.d);

msaisir = saisir_mean(saisir);
msaisir.d = msaisir.d-msaisir.d(1);

% coeff = sum(msaisir.d.^2,2)/K;
[~,i1]=min(abs(str2num(msaisir.v)-3100.0));
[~,i2]=min(abs(str2num(msaisir.v)-900.0));
msaisirsel=selectcol(msaisir,[i1:i2]);
[coeff,ind] = max(msaisirsel.d);
peak = str2num(msaisirsel.v(ind,:));

refcorr.d = msaisir.d./coeff;
refcorr.v = msaisir.v;
refcorr.i = 'reference';

% figure
% plot(str2num(msaisir.v),msaisir.d);
% set(gca,'XDir','reverse');
% set(gcf,'Color',[1 1 1]);
% axis tight;
% figure
% plot(str2num(refcorr.v),refcorr.d);
% set(gca,'XDir','reverse');
% set(gcf,'Color',[1 1 1]);
% axis tight;
end