function [] = Save_Data2GUI(Coeffs, Dict, sesInfo )
%   This fucntion take the data and save them in a format legible by the
% GIFT tool display automatilcally
%   14 Mar 2016                                              M. Morante
%--------------------------------------------------------------------------
%       INPUT:
%   Coeffs [numComp x Mask_Space] -> The spatial maps of each component
%   Dict   [numTime x numComp]    -> The time courses of each component
%   sesInfo -> The main info from the parameters file
%--------------------------------------------------------------------------
%       OUTPUT:
%   Pref_ica_c1_1.mat                   file
%   Pref_ica_br1.mat                    file
%   Pref_sub01_component_ica_s1_.mat    file
%   Pref_sub01_component_ica_s1_.nii    file
%   Pref_sub01_timecourses_ica_s1_.nii  file
%--------------------------------------------------------------------------

    % Main parameters
    Pref = sesInfo.userInput.prefix;
    NoC  = sesInfo.numComp;
    NoT  = sesInfo.diffTimePoints(1);
    
    outPath = [sesInfo.outputDir filesep];
    Mask    = sesInfo.userInput.mask_ind;
    
    
    %%% SAVE Pref_ica_c1_1.mat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Pareameters
    ic = Coeffs;
    tc = Dict';
    
    save([outPath Pref '_ica_c1-1.mat'], 'ic', 'tc');
    
    clear('ic', 'tc');
    
    %%% SAVE Pref_ica_br1.mat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parameters
    compSet.ic = Coeffs;
    compSet.tc= Dict';
    
    save([outPath Pref '_ica_br1.mat'], 'compSet');
    
    clear('compSet');
    
    %%% SAVE Pref_sub01_component_ica_s1_.mat %%%%%%%%%%%%%%%%%%%%%%
    
    % Parameters
    mat = sesInfo.HInfo.V.mat;
    
    save([outPath Pref '_sub01_component_ica_s1_'], 'mat');
    
    clear('mat');
    
    %%% SAVE NIFTI FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %-----% SPATIAL COMPONENTS %-------------------------------------
    
    % Open the .nii parameters template from the mask
    nii = load_nii([outPath Pref 'Mask.hdr']);
    
    % Change the file type
    nii.filetype = 2;
    
    % Actualize new dimmensions
    dim    = ones(1,8);
    dim(1) = 4;
    for i=2:4
        dim(i) = sesInfo.HInfo.DIM(i-1);
    end
    dim(5) = NoC;
    
    nii.hdr.dime.dim = dim;
    
    % Spatial dimensions
    dimX = dim(2);
    dimY = dim(3);
    dimZ = dim(4);
    
    clear('dim');
    
    % Put data
    img = zeros(dimX,dimY,dimZ,NoC);
    
    for i=1:NoC
        aux = zeros(1,dimX*dimY*dimZ);
        if size(Mask,1)==size(Coeffs,2)
            aux(Mask)= Coeffs(i,:);
        else
            aux= Coeffs(i,:); 
        end        
        img(:,:,:,i) = reshape(aux,dimX,dimY,dimZ);
    end
    
%     for i=4:13
%         temp=sprintf('C:\\Users\\hatzi\\Desktop\\test2\\cope3_%d.nii.gz',i);
%         tempii=load_untouch_nii(temp);
%         img(:,:,:,i-3)=double(tempii.img);
%     end
% %     
    for i=4:14
        temp=sprintf('C:\\Users\\hatzi\\Desktop\\test1\\cope%d.nii.gz',i);
        tempii=load_untouch_nii(temp);
        img(:,:,:,i-3)=double(tempii.img);
    end
%     
%     temp=sprintf('C:\\Users\\hatzi\\Desktop\\cope13.nii.gz');
%     tempii=load_untouch_nii(temp);
%     img(:,:,:,1)=double(tempii.img);
    clear('aux', 'dimX', 'dimY', 'dimZ');
   
    
    % Actualize the new img in the nii file
    nii.img = img;
    
    clear('img');
    
    % Actualize the new name
    name = [Pref '_sub01_component_ica_s1_'];
    nii.fileprefix = name;  
    nii.ext = '';
    
    % Save
    path = [outPath name '.nii'];
    save_nii(nii,path);
    
    clear('name', 'path');
    clear('nii');
    
    %-----% TEMPORAL COMPONENTS %-------------------------------------
    
    % Opent the .nii parameters template from the mask
    nii = load_nii([outPath Pref 'Mask.hdr']);
    
    % Change file type
    nii.filetype = 2;
    nii.ext = '';
    
    % Actualize new dimensions
    dim    = ones(1,8);
    dim(1) = 2;
    dim(2) = NoT;
    dim(3) = NoC;
    
    nii.hdr.dime.dim = dim;
    
    clear('dim');
    
    % Actualize data 
    nii.img = Dict;
    
    % Actualize the name of the file
    name = [Pref '_sub01_timecourses_ica_s1_'];
    nii.fileprefix = name;
    
    path = [outPath name '.nii'];
    save_nii(nii,path);
    
    clear('name', 'path');
    clear('nii');
    
    display('SUCCESSFULLY SAVED!          \(^.^)/');
    
end

