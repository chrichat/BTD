function [icasig,Bnew,Cnew] = tpica(T,F,Iter,varargin)

%Tensor PICA
%
% The Iterative TPICA algorithm of Beckman & Smith 2005 for extracting
% F source signals from 3-way fMRI data T of the form I voxels x J time
% points x K subjects. R is the order of decomposition (Number of independent components to
% be estimated). In the initialiazation Yk is the k-th's subject I
% voxels x Jtime points data matrix which is assumed to be PREWHITENED
% (have zero mean columns and white Gaussian noise). FastICA is beeing used
%
% You have to define if the algorithm to be used will be iterative or not
% 
%Iter specifies if the iteartive TPICa will be used 1 or not 0
%


[voxels,timep,subj] = size(T);
Xp=reshape(T,voxels,timep*subj); %reshape 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the mean and check the data
Xpn=Xp;
% [Xpn, mixedmean] = remmean(Xp);
%%%%%%%%%%

[U,D,V] = svd(Xpn);
sn=(sum(diag(D.^2))-sum(diag(D(1:F,:).^2)))/(voxels*(timep*subj-F)); %compute noise variance
W1=V(:,1:F);
W2=D(1:F,1:F)^2/voxels-sn*eye(F);
Xt=Xpn*V(:,1:F)*W2^(-1/2);
Aold=rand(voxels,F);
Bold=rand(timep,F);
Cold=rand(subj,F);
e=10^-7;
Bnew=Bold;
Cnew=Cold;
Rold=rand(F,F);
flag=1;

while flag==1
%     [icasig,Rnew, Wnn] = fastica(Xt','initGuess',Rold);     
    [icasig,Rnew, Wnn] = icaML(Xt',F);
    Anew=icasig';
    M2=W1*W2^0.5*(Rnew/(Rnew'*Rnew)); %% better approximation of Rnew theoretically Rnew'*Rnew is the identity
    for i=1:size(icasig,1)
        zf=reshape(M2(:,i),timep,subj);
        [uf,df,vf] = svd(zf);
        Bnew(:,i)=uf(:,1);
        Cnew(:,i)=df(1,1)*vf(:,1);    
    end
    if Iter==0
        flag=0;
    else
        delta=norm(Anew-Aold)+norm(Bnew-Bold)+norm(Cnew-Cold);
        if delta<e
            flag=0;
        else
           Cold=Cnew;
           Bold=Bnew;
           Rold=W2^-0.5*V(:,1:F)'*(kr(Cnew,Bnew));
        end
    end
end



 




