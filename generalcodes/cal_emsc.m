function [ZCorrected,ZResiduals,ZParameters]=cal_emsc(ZSaisir,EMSCModel,ZWeights,option);
%cal_emsc  [ZCorrected,ZResiduals,ZParameters]=cal_emsc(ZSaisir,EMSCModel,ZWeights,option)
%
%           'runs EMSC on ZSaisir using the EMSCModel defined by make_emsc_modfunc
%           model functions can be added to EMSCModel by add_spec_to_mod'
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
%  First version: 12.03.05                                                             %
%  Revision:                                                                           %
%              One revision for speed improvement by Kristian Liland                   %
%              28.10.2016 option added for standard EMSC correction and                %
%                         correction of only additive effects                          %
%                         complete revision
%                                                                                      %
%                                                                                      %
%                                                                                      %
%--------------------------------------------------------------------------------------%
%  function [ZCorrected,ZResiduals,ZParameters]=run_emsc(ZSaisir,EMSCModel,ZWeights)   %
%                                                                                      %
%  Runs emsc using emsc basic model fucntions and constituent difference spectra       %
%                                                                                      %
%                                                                                      %
%  Input:   ZSaisir: saisire structure                                                 %
%           EMSC modell functions (in EMSC structure)                                  %                      
%           ZWeights: saisire structure                                                %
%  option=0: only correction of additive terms                                         %
%                                                                                      %
%  If option is not defined or 1: default is EMSC with correction of additive and      %
%                                 multiplicative variations                            %   
%                                                                                      %
%                                                                                      %
%                                                                                      %
%  Output:  saisir structures for ZCorrected,ZResiduals,ZParameters                    %
%                                                                                      %
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isempty(option)) 
    Option=1;    
else
    Option=option;
end

%% check if model functions and spectra have the same dimensions
    [NZ, MZ]=size(ZSaisir.d);
    [NE, ~]=size(EMSCModel.ModelVariables);
    if (MZ~=NE)
        error('model functions and spectra do not have the same dimensions')
    end

%% use weights in case they are defined:
    if (isempty(ZWeights)) 
        ZSaisirSel=ZSaisir;        
    else
        ZSaisirSel=ZSaisir;
        [WN,~]=size(ZWeights.d);
        if (WN>1)
            for i=1:NZ
                ZSaisirSel.d(i,:)=ZSaisir.d(i,:).*ZWeights.d(2,:);
            end
        else
            for i=1:NZ
                ZSaisirSel.d(i,:)=ZSaisir.d(i,:).*ZWeights.d;
            end
        end
    end    
    
%% initialise and weight the model    
    [Nx Ny]=size(ZSaisir.d);
    Residuals=zeros(Nx,Ny);
    Corrected=zeros(Nx,Ny);
    [~, M]=size(EMSCModel.Model);
    Parameters=zeros(Nx,M);
    Model=EMSCModel.Model;
    if (isempty(ZWeights)) 
        Model=EMSCModel.Model;
    else
        [WN,~]=size(ZWeights.d);
        if (WN>1)
            for i=1:M
                Model(:,i)=EMSCModel.Model(:,i).*ZWeights.d(i,:)';
            end
        else
            for i=1:M
                Model(:,i)=EMSCModel.Model(:,i).*ZWeights.d';
            end
        end
    end   
    
    
% Handle all spectra in one go to save time
PP = all(ZSaisirSel.d==0,2);  % Finds all spectra which are zero

if any(PP) % if there are any spectra that are zero, then we treat only the ones that are not zero, i.e. ~PP
    Parameters0  = (Model\ZSaisirSel.d(~PP,:)')';    % Test this with respect to least squares
    
     % Remove baseline
    k=(1:(EMSCModel.NumBasicModelFunc-1));
    Corrected0 = ZSaisir.d(~PP,:)-Parameters0(:,k)*EMSCModel.Model(:,k)';
    if EMSCModel.NumBasicModelFunc == 1
        Corrected0 = ZSaisir.d-Parameters0*EMSCModel.Model';
    end
    % Remove bad spectra
    if (EMSCModel.NbadSpec>0)
        k=(EMSCModel.NumBasicModelFunc+EMSCModel.NgoodSpec+1):M;   % bad spectra
        Corrected0 = Corrected0-Parameters0(:,k)*EMSCModel.Model(:,k)';
    end
    
    % Correct multiplicative effect
    if (Option==1)
        Corrected0=bsxfun(@times,Corrected0,1./Parameters0(:,EMSCModel.NumBasicModelFunc)); % correct multipl. eff.
    end
    
    % Calculate the residuals
    k=(1:M);
    Residuals0=ZSaisir.d(~PP,:)-Parameters0(:,k)*EMSCModel.Model(:,k)';
    
    Corrected(PP,:)   = ZSaisir.d(PP,:);
    Corrected(~PP,:)  = Corrected0;
    Residuals(~PP,:)  = Residuals0;
    Residuals(PP,:)   = 0;
    Parameters(~PP,:) = Parameters0;
    Parameters(PP,:)  = 0;
    
    ZCorrected.d=Corrected;
    ZCorrected.v=ZSaisir.v;
    ZCorrected.i=ZSaisir.i;
    
    ZResiduals.d=Residuals;
    ZResiduals.v=ZSaisir.v;
    ZResiduals.i=ZSaisir.i;
    
    ZParameters.d=Parameters;
    ZParameters.v=EMSCModel.ModelSpecNames;
    ZParameters.i=ZSaisir.i;
else
    Parameters0  = (Model\ZSaisirSel.d')';   
    
    % Remove baseline
    k=1:(EMSCModel.NumBasicModelFunc-1);
    Corrected0 = ZSaisir.d-Parameters0(:,k)*EMSCModel.Model(:,k)';    
    if EMSCModel.NumBasicModelFunc == 1
        Corrected0 = ZSaisir.d-Parameters0*EMSCModel.Model';
    end
    
    % Remove bad spectra
    if (EMSCModel.NbadSpec>0)
        k=(EMSCModel.NumBasicModelFunc+EMSCModel.NgoodSpec+1):M;  % bad spectra
        Corrected0 = Corrected0-Parameters0(:,k)*EMSCModel.Model(:,k)';
    end
    
    % Correct multiplicative effect
    if (Option==1)
        Corrected0=bsxfun(@times,Corrected0,1./Parameters0(:,EMSCModel.NumBasicModelFunc));
    end
    
    % Calculate the residuals
    k=1:M;
    Residuals0=ZSaisir.d-Parameters0(:,k)*EMSCModel.Model(:,k)';
    
    ZCorrected.d=Corrected0;
    ZCorrected.v=ZSaisir.v;
    ZCorrected.i=ZSaisir.i;
    
    ZResiduals.d=Residuals0;
    ZResiduals.v=ZSaisir.v;
    ZResiduals.i=ZSaisir.i;
    
    ZParameters.d=Parameters0;
    ZParameters.v=EMSCModel.ModelSpecNames;
    ZParameters.i=ZSaisir.i;
end

end