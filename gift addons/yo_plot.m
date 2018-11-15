%--------------------------------------------------------------------------
% This function complete the whole porces to extract and to save the data
%  in the middle the data are analyzed
%--------------------------------------------------------------------------   
%       INTPUS
%   Pref -> The prefix of the parameters
%   Dir  -> The path to the parameter file
%--------------------------------------------------------------------------
function [Coeffs, Dict, W]=yo_plot(name)
%     % Adding pahts
%     Pref_ind=strfind(name,'_'); %identify the length of the prefix based on the presence of _ 
%     Pref=name(1:Pref_ind-1);
%     
%     addpath(genpath(['matlabtoolsforplottingfmri' filesep]))
%     addpath(genpath(['FastICA' filesep]))
%     addpath(genpath(['Auxiliars' filesep]))



%% READ PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     name = [Dir filesep Pref '_ica_parameter_info.mat'];
    load(name);
    
    % Is initialiced
    if (sesInfo.isInitialized == 1)
        display ('WARNING: The parametes have been analyzed before...      ¯\(O_O)')
    else
        sesInfo = Mod_Param(sesInfo); % Modification of the parameters
    end
    
    %Extract important information
    Pref = sesInfo.userInput.prefix;
    NoC  = sesInfo.numComp;
    NoT  = sesInfo.diffTimePoints;
    
    path_name_nii = sesInfo.HInfo.V.fname;
    nii = load_untouch_nii(path_name_nii); 
    
    % Exporting 4d image 
    Y = nii.img; 
    
    clear('nii');
    
    
%%  RUN ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% (Yannis part...)
    preProcType='remove mean per timepoint';
    Y = icatb_preproc_data(Y, preProcType);
    %Build the data matrix applaying the mask
    TIM = sesInfo.userInput.diffTimePoints;
    for n = 1:TIM
        
        % UNFOLD FUNCTION: For each time n
        Vo(n,:) = M_UnFold(Y(:,:,:,n),sesInfo.userInput.mask_ind);

    end
    
    % FOLD FUNCTION
    fold = M_Fold(Vo(n,:),sesInfo.HInfo.DIM,sesInfo.userInput.mask_ind);
%    NoC=20;
    % Applying ICA % You only have to change our method here
    nICs = NoC; % <======================================================== I put the components from the parameters file
   % [icasig1, Dict, W] = fastica(X, 'verbose', 'on', 'numOfIC', nICs);
    
    nPCs = NoC; % <======================================================== Idem here !
    rng(211,'twister');
    [Coeffs, Dict, W] = fastica(Vo, 'verbose', 'on', 'numOfIC', nICs, 'firstEig', 1, 'lastEig', nPCs); 


% [CCC,CP]=corr(design_mat,Dict);      % calculate correlations and p-values
% [CCCm,AAAm]=max(abs(CCC)');
for cmp = 1:size(Dict,2)
   varcmp = var(Coeffs(cmp,:));
   Coeffs(cmp,:) = Coeffs(cmp,:)/sqrt(varcmp);
   signcmp = sign(max(Coeffs(cmp,:)));
   Coeffs(cmp,:) = signcmp*Coeffs(cmp,:);
   
   Dict(:,cmp) = signcmp*Dict(:,cmp)*sqrt(varcmp);
    
end

%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Save_Data2GUI(Coeffs,Dict,sesInfo);
   

disp('');