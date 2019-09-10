function CVPar = CVParGenerator(CVPar,X)
% function CV = CVParGenerator(NX,CVType,CVFold)
% Fill in the dscription!!!
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