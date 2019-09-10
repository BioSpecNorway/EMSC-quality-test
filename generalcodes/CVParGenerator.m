function CVPar = CVParGenerator(CVPar,X)
% function CV = CVParGenerator(CVPar,X)
% Input: CVPar contains CVPar.CVType which splits the data by
%        'Full CV' or leave-one-out
%        'Contiguous Blocks' in K blocks, number specified in CVPar.cvfold
%        'Venetian Blinds' in K blocks, number specified in CVPar.cvfold
%        'Random Subsets' in K blocks, number specified in CVPar.cvfold
%        'Specified' according to the group identified by the labels in
%                    name position in CVPar.cvind (1x2)
% Output: CVPar: cv containing labels of groups
%                fold number of blocks
%                ind index in name for *Specified*
%                
% by V.Tafintseva 20.05.2015

%%
CVType = CVPar.CVType;

if strcmp(CVType,'Specified')
    cvind = CVPar.cvind;
    Gr = create_group1(X,cvind(1),cvind(end));
    CVPar.cvfold = max(Gr.d);
    for i=1:max(Gr.d)
    IndexOut = seekstring(X.i(:,cvind(1):cvind(end)),Gr.g.i(i,:));
    CV(IndexOut,1) = i;
    end
else
CVFold = CVPar.cvfold;
NX = size(X.d,1);

CV = [];
NSamples = round(NX/CVFold);

if strcmp(CVType,'Full CV') 
    CV = (1:NX)';
    CVPar.cvind = []; CVPar.cvfold = NX;
elseif strcmp(CVType,'Contiguous Blocks')
    for i = 1:CVFold-1
    CV = [CV; repmat(i,NSamples,1)];
    CVPar.cvind = [];
    end
    CV = [CV; repmat(CVFold,NX-i*NSamples,1)];
    CVPar.cvind = [];
elseif strcmp(CVType,'Random Subsets')
    CV = randi(CVFold,NX,1);    
elseif strcmp(CVType,'Venetian Blinds')
    CV = [CV; repmat((1:CVFold)',NSamples-1,1)];
    CVPar.cvind = [];
    if NX-(NSamples-1)*CVFold > CVFold
        CV = [CV; (1:CVFold)'; (1:NX-(NSamples)*CVFold)'];
    else 
        CV = [CV; (1:NX-(NSamples-1)*CVFold)'];
    end
else error('Unknown CV type')
end % if

end
CVPar.cv = CV;

end