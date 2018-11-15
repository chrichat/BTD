    
function  plotting_tpica(d,mask,t)

    T=permute(d,[2,1,3]);
    F=4;
    
    [voxels,timep,subj] = size(T);
    Xp=reshape(T,voxels,timep*subj); %reshape
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove the mean and check the data
    
    [Xpn, mixedmean] = remmean(Xp);
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
    [icasig,Rnew, W] = fastica(Xt', 'maxNumIterations' ,10000,'initGuess',Rold);
    
    Anew=icasig';
    M2=W1*W2^0.5*(Rnew/(Rnew'*Rnew)); %% better approximation of Rnew theoretically Rnew'*Rnew is the identity
    for i=1:F
        zf=reshape(M2(:,i),timep,subj);
        [uf,df,vf] = svd(zf);
        Bnew(:,i)=uf(:,1);
        Cnew(:,i)=df(1,1)*vf(:,1);    
    end
    
    mask1=mask(1,1:4096)';

    z11=icasig(1,1:4096);
    z12=icasig(1,4097:8192);
    z13=icasig(1,8193:12288);
    z11(~mask1)=0;
    z12(~mask1)=0;
    z13(~mask1)=0;
    
    z21=icasig(2,1:4096);
    z22=icasig(2,4097:8192);
    z23=icasig(2,8193:12288);
    z21(~mask1)=0;
    z22(~mask1)=0;
    z23(~mask1)=0;

    z31=icasig(3,1:4096);
    z32=icasig(3,4097:8192);
    z33=icasig(3,8193:12288);
    z31(~mask1)=0;
    z32(~mask1)=0;
    z33(~mask1)=0;
    
    z41=icasig(4,1:4096);
    z42=icasig(4,4097:8192);
    z43=icasig(4,8193:12288);
    z41(~mask1)=0;
    z42(~mask1)=0;
    z43(~mask1)=0;
    
    bottom1 = min([min(min(z11)),min(min(z12)),min(min(z13)),min(min(z21)),min(min(z22)),min(min(z23)),min(min(z31)),min(min(z32)),min(min(z33))]);
    top1  = max([max(max(z11)),max(max(z12)),max(max(z13)),max(max(z21)),max(max(z22)),max(max(z23)),max(max(z31)),max(max(z32)),max(max(z33))]);

    if abs(bottom1)>abs(top1)
        top1=abs(bottom1);
    end
    mask2=0.4*top1*mask1';
    mask2(~mask1)=0;

    figure
    subplot(4,5,1)
    imagesc(imrotate(reshape(abs(z11)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);
    subplot(4,5,2)
    imagesc(imrotate(reshape(abs(z12)+mask2,64,64),90))
    title('Spatial Mode 1','fontweight','bold','fontsize',16);
    caxis manual
    caxis([0 top1]);
    subplot(4,5,3)
    imagesc(imrotate(reshape(abs(z13)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);
    subplot(4,5,[4,5])
    plot(Bnew(:,1));
    ylabel(' % Bold change')
    title('Temporal Mode 1','fontweight','bold','fontsize',10);
    subplot(4,5,6)
    imagesc(imrotate(reshape(abs(z21)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);
    subplot(4,5,7)
    imagesc(imrotate(reshape(abs(z22)+mask2,64,64),90))
    title('Spatial Mode 2','fontweight','bold','fontsize',16);
    caxis manual
    caxis([0 top1]);
    subplot(4,5,8)
    imagesc(imrotate(reshape(abs(z23)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);
    subplot(4,5,[9,10])
    plot(Bnew(:,2));
    ylabel(' % Bold change')
    title('Temporal Mode 2','fontweight','bold','fontsize',10);
    subplot(4,5,11)
    imagesc(imrotate(reshape(abs(z31)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);
    subplot(4,5,12)
    imagesc(imrotate(reshape(abs(z32)+mask2,64,64),90))
    title('Spatial Mode 3','fontweight','bold','fontsize',16);
    caxis manual
    caxis([0 top1]);
    subplot(4,5,13)
    imagesc(imrotate(reshape(abs(z33)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);    
    subplot(4,5,[14,15])
    plot(Bnew(:,3));
    ylabel(' % Bold change')
    title('Temporal Mode 3','fontweight','bold','fontsize',10);
    subplot(4,5,16)
    imagesc(imrotate(reshape(abs(z41)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);
    subplot(4,5,17)
    imagesc(imrotate(reshape(abs(z42)+mask2,64,64),90))
    title('Spatial Mode 3','fontweight','bold','fontsize',16);
    caxis manual
    caxis([0 top1]);
    subplot(4,5,18)
    imagesc(imrotate(reshape(abs(z43)+mask2,64,64),90))
    caxis manual
    caxis([0 top1]);
    subplot(4,5,[19,20])
    plot(Bnew(:,4));
    ylabel(' % Bold change')
    title('Temporal Mode 4','fontweight','bold','fontsize',10);
   

end
