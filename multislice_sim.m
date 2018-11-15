function [corr_tpica, corr_cpd, corr_btd, corr_pfac2, corr_btd2]=multislice_sim(SNR, dataset,iter)

% This function reproduces the results presented in section 3.1.1 of [1]
%
%Input -SNR:level of Signal to Noise Ratio, in the paper have been used
%           values from 0.06 to 0.15. Signal to noise ratio is defined as the
%           frobenius norm of the signal to the frobenius norm of the noise
%     -dataset:string. Which dataset will be used. The options are 'A',
%           'G','C' and 'I'.
%      -iter:Number of iterations for every methond (iter=30 used in paper)
%Output -The mean correlation among sources for every method per iteration


% Part of the code used in the script for the generation of the data has
% been provided by A. Stegeman and was presented in [2].
% In order to run the script tensorlab 3.0 is nedded [3].
% The function lp_nls has been provided by Ot. Debals and its use is
% presented in [4], the use of it fells  under the same license as
% tensorlab 3.0 [3].
%
% [1]C. Chatzichristos et al.'Blind fMRI Source Unmixing via Higher-Order Tensor Decompositions'
% [2]A. Stegeman, 'Comparing Independent Component Analysis and the PARAFAC
% model for artifcial multisubject fMRI data'
% [3]N. Vervliet et al. 'Tensorlab 3.0' Available: http://www.tensorlab.net
% [4]M. Bousse et al.'A tensor-based method for large-scale blind source separation using segmentation'



if dataset=='A'
    gen_data_A %Function provided by A. Stegeman
elseif dataset=='C'
    gen_data_C %Function provided by A. Stegeman
elseif dataset=='G'
    gen_data_G %Function provided by A. Stegeman
elseif dataset=='I'
    gen_data_I
else
    error('Dataset %s is not an option. The available datasets are A,G,C and I\n',dataset);
end

datasetA                  %Function provided by A. Stegeman
Ub{1,1}(:,1)=f(mask==1);  % Masking the in brain voxels
Ub{1,1}(:,2)=s(mask==1);
Ub{1,1}(:,3)=l(mask==1);
if size(B0,2)>size(C0,1)
    Ub{1,2}=mean(reshape(B0,size(B0,1),3,[]),3);    %Dataset C and I
else
    Ub{1,2}=B0;
end


    
for kk=1:iter
  
    %SNR decrease%        
    c=norm(Xpreps,'fro')/(SNR*(norm(Xprepn,'fro')));
    d=reshape(Xpreps+c*Xprepn,size(Xprepn,1),[],3);

    %Remove mean%
    d=permute(d,[2,1,3]);   
    Xpn1=reshape(d,size(d,1),[]);
    [Xpn,menv]=remmean(Xpn1);
    
    %CPD%
    d=reshape(Xpn,size(d));
    R=3;
    tic
    options = struct;
    options.Initialization = @cpd_gevd;
    options.Algorithm = @cpd_nls;
    % options.AlgorithmOptions.LineSearch = @cpd_lsb; % Add exact line search.
    options.Compression = @mlsvd;
    options.ExploitStructure=false;
    options.Refinement=false; % Remove in case you want refinement it gets much slower and without significant gain
    [U,output] = cpd(d,R,[196*3,196,3],options,'Display',1);
    toc
    for i=1:size(Xprep_3d,3)
        Ucpd{1,1}(:,i)=U{1,1}(mask==1,i);
    end
    Ucpd{1,2}=U{1,2};
    [ corr_cpd(kk),scorr_cpd,tcorr_cpd ] = mean_corr( Ucpd,Ub );  % Corr_cpd is the mean cpd used THE MAIN METRIC
    means_cpd(kk)=mean(scorr_cpd);
    meant_cpd(kk)=mean(tcorr_cpd);

    %%%%%%%%%%%TPICA %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    Sica=d;
    Xpn=reshape(Sica,size(Sica,1),[]);
    [voxels,timep,subj] = size(Sica);
    [~,D,V] = svd(Xpn);
    sn=(sum(diag(D.^2))-sum(diag(D(1:R,:).^2)))/(voxels*(timep*subj-R)); %compute noise variance
    W1=V(:,1:R);
    W2=D(1:R,1:R)^2/voxels-sn*eye(R);
    Xt=Xpn*V(:,1:R)*W2^(-1/2);
    Bnew=rand(timep,R);
    Cnew=rand(subj,R);
    Rold=rand(R,R);
    [icasig,Rnew, ~] = fastica(Xt', 'maxNumIterations' ,10000,'initGuess',Rold);
    Anew=icasig';
    M2=W1*W2^0.5*(Rnew/(Rnew'*Rnew)); %% better approximation of Rnew theoretically Rnew'*Rnew is the identity
    for i=1:R
        zf=reshape(M2(:,i),timep,subj);
        [uf,df,vf] = svd(zf);
        Bnew(:,i)=uf(:,1);
        Cnew(:,i)=df(1,1)*vf(:,1);    
    end
    W1=kr(Cnew,Bnew)';
    icasig = W1 * Xpn1';
    for i=1:size(Xprep_3d,3)
        Utpica{1,1}(:,i)=icasig(i,mask==1);
    end
    Utpica{1,2}=Bnew;
    Utpica{1,3}=Cnew;
%     Utpica=cellfun(@(u,v)v*u,Utpica,Vica,'UniformOutput',false); %%uncompressing results
    Utpica=Utpica(1,1:2);
    toc
    
    [corr_tpica(kk),scorr_tpica,tcorr_tpica] = mean_corr( Utpica,Ub ); % Corr_tpica is the mean cpd used THE MAIN METRIC
    means_tpica(kk)=mean(scorr_tpica);
    meant_tpica(kk)=mean(tcorr_tpica);
    
    %%%%%%%%%%% BTD %%%%%%%%%%%%% 
    L=10;   %Change L accrdingly
    N = [2 2]; % Order of the tensors
    T4d=reshape(d,64,64*3,size(d,2),[]);
    U=cpd_gevd(d,R);
    U{1,1}=icasig';
    U{1,2}=Utpica{1,2};
    clear U_gevd
     for k=1:R
       tempimg=reshape(U{1,1}(:,k),64,[]);
       [U1,S,V] = svd(tempimg);
       U_gevd{1,k}{1,1}=U1(:,1:L)*S(1:L,1:L);
       U_gevd{1,k}{1,2}=V(:,1:L);
       U_gevd{1,k}{1,4}=U{1,3}(:,k);
       U_gevd{1,k}{1,3}=U{1,2}(:,k);
     end
    tic
    [V,S,~] = mlsvd(T4d,[30,30,196,3]); %compression
    for k=1:size(U_gevd,2)
        U2{1,k}= cellfun(@(u,v)v'*u,U_gevd{1,k}(1,1:4),V,'UniformOutput',false); %%compressing the initializations with the same compression matrix V
        U2{1,k}{1,5}=diag(ones(1,L));
    end
    U_comp = lp_nls(S,U2,N); %Provided by Otto debals, 
    toc
    for k=1:size(U_gevd,2)
        U21{1,k}=cellfun(@(u,v)v*u,U_comp{1,k}(1,1:4),V,'UniformOutput',false); %%uncompressing results
        tempimg=U21{1,k}{1,1}*U21{1,k}{1,2}';
        tempimg=reshape(tempimg,1,[]);
        Ubtd{1,1}(:,k)=tempimg(mask==1);
        Ubtd{1,2}(:,k)=U21{1,k}{1,3};
    end
    [corr_btd(kk),scorr_btd,tcorr_btd] = mean_corr( Ubtd,Ub );
    means_btd(kk)=mean(scorr_btd);
    meant_btd(kk)=mean(tcorr_btd);
    
    [Vp2,Sp2,~] = mlsvd(d,[196*3,196,3]); %compression
    [A,H,C,P,~,~]=parafac2(d,3);
    Upfac2{1,1}=A;
    Upfac2{1,2}=mean(reshape([P{1,1}*H,P{1,2}*H,P{1,3}*H],size(B0,1),3,[]),3);
%     Upfac2{1,2}=P{1,3}*H;
    Upfac2{1,3}=C;
    Upfac2=cellfun(@(u,v)v*u,Upfac2,Vp2,'UniformOutput',false); %%uncompressing results
    for k=1:3
        temp=Upfac2{1,1}(:,3);
        Upfac22{1,1}(:,k)=temp(mask==1);
    end
    Upfac22{1,2}=Upfac2{1,2};
%     Upfac22{1,2}=icatb_filt_data(Upfac22{1,2},1,0.8);
    Upfac22{1,3}=Upfac2{1,3};
    [corr_pfac2(kk),scorr_pfac2,tcorr_pfac2] = mean_corr( Upfac22,Ub );
    means_pfac2(kk)=mean(scorr_pfac2);
    meant_pfac2(kk)=mean(tcorr_pfac2);
    
    
    %%%%%%%%%%%%BTD2-NLS%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K=size(d,3);
    H = eye(R);    
    for k = 1:K
        Qk= data(:,:,k)'*(A*diag(C(k,:))*H');
        P{k}= Qk*psqrt(Qk'*Qk);
        T4(:,:,k) = data(:,:,k)*P{k};
    end
    iter2=50;    
    for j=1:iter2
        U_btd2 = lp_nls(T4,U_btd2,N);
        if i<R
            for i=1:R
                aa(:,:,i)=U_btd2{1,i}{1,1}*U_btd2{1,i}{1,2}';
                C(:,i)=U_btd2{1,i}{1,3};
                H(:,i)=U_btd2{1,i}{1,4};
            end
            A=reshape(A,[],3);
            for k = 1:K
                Qk= data(:,:,k)'*(A*diag(C(k,:))*H');
                P{k}= Qk*psqrt(Qk'*Qk);
                T44(:,:,k) = data(:,:,k)*P{k};
            end
            T4=reshape(T44,size(T4d,1),size(T4d,2),size(T4,3),[]);
        end
        endp
    
    clear A
    for i=1:R
        A(:,i)=squeeze(reshape(U_btd2{1,i}{1,1}*U_btd2{1,i}{1,2}',[],1));
    end
    plotme(A',size(ymask,2),size(xmask,2))

    
end
end