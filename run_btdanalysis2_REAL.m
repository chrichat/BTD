%Script for running CPD,TPICA,BTD, PARAFAC2 and BTD2



% a=load_nii('H:\Human Connectome\101915\MNINonLinear\Results\tfMRI_RELATIONAL_RL\tfMRI_RELATIONAL_RL.nii.gz');
% a.filetype=1;
% save_nii(a, 'H:\Connectome- Relational\101915.nii')
% a=load_nii('H:\Human Connectome\105115\MNINonLinear\Results\tfMRI_RELATIONAL_RL\tfMRI_RELATIONAL_RL.nii.gz');
% a.filetype=1;
% save_nii(a, 'H:\Connectome- Relational\105115.nii')
% 
% a=load_nii('H:\Human Connectome\103111\MNINonLinear\Results\tfMRI_RELATIONAL_RL\tfMRI_RELATIONAL_RL.nii.gz');
% a.filetype=1;
% save_nii(a, 'H:\Connectome- Relational\100311.nii')
% 
% a=load_nii('H:\Human Connectome\133928\MNINonLinear\Results\tfMRI_RELATIONAL_RL\tfMRI_RELATIONAL_RL.nii.gz');
% a.filetype=1;
% save_nii(a, 'H:\Connectome- Relational\110411.nii')
% 
% a=load_nii('H:\Human Connectome\117122\MNINonLinear\Results\tfMRI_RELATIONAL_RL\tfMRI_RELATIONAL_RL.nii.gz');
% a.filetype=1;
% save_nii(a, 'H:\Connectome- Relational\117122.nii')

j=[5,6,7,11,12,13,14,15,17,19,20,21,22,23,27,28,29,30];
for i=1:size(j,2)
    if j(i)<10
        temp=sprintf('U:\\OpenfMRI-ds157\\ds157\\sub-0%d\\func\\bold_fw8.feat\\filtered_func_data.nii.gz',j(i));
    else
        temp=sprintf('U:\\OpenfMRI-ds157\\ds157\\sub-%d\\func\\bold_fw8.feat\\filtered_func_data.nii.gz',j(i));
    end
    d(i)=load_nii(temp);
end



%%%%
for i=1:size(d,2)
%     dat(:,:,:,:,i)=d(i).img(minx:maxx,miny:maxy,minz:maxz,:);
    dat(:,:,:,:,i)=double(d(i).img);
end
    
datas=double(reshape(dat,[],size(dat,4),size(dat,5)));
% Es=randn(size(datas));
% datas=datas+0.01*max(max(max(datas)))*Es;
for i=1:size(d,2)
%     if i<7
%         datas(:,:,i)=datas(:,:,i)+comp1;    %Augment data
%     elseif i>11
%         datas(:,:,i)=datas(:,:,i)+comp3;
%     else
%         datas(:,:,i)=datas(:,:,i)+comp2;        
%     end
    [datall(:,:,i),mv]=remmean(datas(:,:,i));
end
R=8;
% datall=permute(datall,[2,1,3]);
% [A,H,C,P,fit,AddiOutput]=parafac2(datall,20);


%%Compression
tic
[Vica,Sica,sv] = mlsvd(datall,[size(datas,2)*size(datas,3),size(datas,2),size(datas,3)]); %%% 3d compression
toc
databtd=reshape(datall,size(dat,1),[],size(datall,2),size(datall,3));
databtd5=reshape(datall,size(dat,1),size(dat,2),[],size(datall,2),size(datall,3));
clear dat
sizs=sort(size(databtd)); %% sort the dimensions use the 2nd smaller for compression, smaller is patients
tic
[Vbtd,Sbtd,svbtd] = mlsvd(databtd,[sizs(2)-20,sizs(3),sizs(3),sizs(1)]); %%% 4d compression for btd
toc


%%%%%%%%%%%TPICA %%%%%%%%%%%%
% Sica=da1;
load('U:\OpenfMRI-ds157\FW_8\test1\cpdtes_ica_parameter_info.mat')
Mask    = sesInfo.userInput.mask_ind;
% dat1=datall(Mask,:,:);

dat1=Sica;
Xpn=reshape(dat1,size(dat1,1),[]);
Xpn=remmean(Xpn);
[voxels,timep,subj] = size(dat1);
[U,D,V] = svd(Xpn,0);
sn=(sum(diag(D.^2))-sum(diag(D(1:R,:).^2)))/(voxels*(timep*subj-R)); %compute noise variance
W1=V(:,1:R);
W2=D(1:R,1:R)^2/voxels-sn*eye(R);
Xt=Xpn*V(:,1:R)*W2^(-1/2);
Aold=rand(voxels,R);
Bold=rand(timep,R);
Cold=rand(subj,R);
e=10^-7;
Bnew=Bold;
Cnew=Cold;
Rold=rand(R,R);
[icasig,Rnew, W] = fastica(Xt', 'maxNumIterations' ,10000,'initGuess',Rold);
Anew=icasig';
M2=W1*W2^0.5*(Rnew/(Rnew'*Rnew)); %% better approximation of Rnew theoretically Rnew'*Rnew is the identity
for i=1:R
    zf=reshape(M2(:,i),timep,subj);
    [uf,df,vf] = svd(zf);
    Bnew(:,i)=uf(:,1);
    Cnew(:,i)=df(1,1)*vf(:,1);    
end
Utpica{1,1}=icasig';
Utpica{1,2}=Bnew;
Utpica{1,3}=Cnew;
Utpica=cellfun(@(u,v)v*u,Utpica,Vica,'UniformOutput',false); %%uncompressing results
y11n=Utpica{1,2};
y11n=icatb_convertToZScores(y11n');
y11n=icatb_filt_data(y11n',1,0.15);
plotme(y11n');
temp=sprintf('U:\\OpenfMRI-ds157\\nofw\\TPICA\\tpica_o_noisy');
save(temp,'Utpica');

ri=1;
for R=5:14
    clear options
    options = struct;
    options.Initialization = @cpd_gevd;
    options.Algorithm = @cpd_als;
    options.Compression = false; %use the compression outside to gain time
    options.ExploitStructure=false;
    options.Refinement=false;
    [Uals1,output] = cpd(Sica,R,[274*5,274,5],options,'Display',1); %use the compression outside to gain time
    Uals=cellfun(@(u,v)v*u,Uals1,Vica,'UniformOutput',false); %%uncompressing results
    y1n=Uals{1,2};
    y1n=icatb_convertToZScores(y1n');
    y1n=icatb_filt_data(y1n',1,0.15);
    temp=sprintf('U:\\OpenfMRI-ds157\\nofw\\CPD\\cpd_o_noisy%d',R);
    save(temp,'Uals');
    plotme(y1n')
    title(['Als',num2str(R)])
    ri=ri+1
end

%%% Find the multilinear Singular Values
[Uen,Sen,sven] = mlsvd(databtd);
figure;
for n = 1:size(Uen,2)
    subplot(1,size(Uen,2),n);
    semilogy(sven{n},'.-');
    xlim([1 length(sven{n})])
end

clearvars -except Sbtd Vbtd

n=1;
R=[8,10,15,20];
L=[20,25,30];
for j=1:size(L,2)
    for i=1:size(R,2)
        U=cpd_gevd(reshape(Sbtd,[],size(Sbtd,3),size(Sbtd,4)),R(i));
             for k=1:R(i)
               tempimg=reshape(U{1,1}(:,k),size(Sbtd,1),[]);
               [U1,S,V] = svd(tempimg);
               U_gevd{1,k}{1,1}=U1(:,1:L(j))*S(1:L(j),1:L(j));
               U_gevd{1,k}{1,2}=V(:,1:L(j));
               U_gevd{1,k}{1,4}=U{1,3}(:,k);
               U_gevd{1,k}{1,3}=U{1,2}(:,k);
               U_gevd{1,k}{1,5}=diag(ones(1,L(j)));
             end
        Ubtd= lp_nls(Sbtd,U_gevd,[2 2]);
        for ik=1:R
            temp=cellfun(@(u,v)v*u,Ubtd{1,ik}(1,1:4),Vbtd,'UniformOutput',false); %%uncompressing results
            Ubtd2{1,1}(:,ik)=reshape(temp{1,1}*temp{1,2}',[],1);
            Ubtd2{1,2}(:,ik)=temp{1,3};
            Ubtd2{1,3}(:,ik)=temp{1,4};
        end
        clear Ubtd
        clear temp
        yn=Ubtd2{1,2};
        yn=icatb_convertToZScores(yn');
        yn=icatb_filt_data(yn',1,0.15);
        temp=sprintf('U:\\OpenfMRI-ds157\\FW_8\\BTD\\btd_o_noisy%d_L%d',R(i),L(j));
        save(temp,'Ubtd2');
        plotme(yn')
        title(['btd',num2str(R(i))])
        clear Ubtd
        n=n+1
    end
end

%%%% 5d%%%
R=8;
L=25;
U=cpd_gevd(reshape(databtd5,[],size(Sbtd,3),size(Sbtd,4)),R);
for k=1:R
    tempimg=reshape(U{1,1}(:,k),size(databtd5,1),size(databtd5,2),[]);
    [U2,S2,V2] = mlsvd(tempimg);
    U_gevd{1,k}{1,1}=U1(:,1:L)*S(1:L,1:L);
    U_gevd{1,k}{1,2}=V(:,1:L);
    U_gevd{1,k}{1,4}=U{1,3}(:,k);
    U_gevd{1,k}{1,3}=U{1,2}(:,k);
    U_gevd{1,k}{1,5}=diag(ones(1,L));
end
Ubtd= lp_nls(Sbtd,U_gevd,[2 2]);
for ik=1:R
    temp=cellfun(@(u,v)v*u,Ubtd{1,ik}(1,1:4),Vbtd,'UniformOutput',false); %%uncompressing results
    Ubtd2{1,1}(:,ik)=reshape(temp{1,1}*temp{1,2}',[],1);
    Ubtd2{1,2}(:,ik)=temp{1,3};
    Ubtd2{1,3}(:,ik)=temp{1,4};
end
clear Ubtd
clear temp
yn=Ubtd2{1,2};
yn=icatb_convertToZScores(yn');
yn=icatb_filt_data(yn',1,0.15);
temp=sprintf('U:\\OpenfMRI-ds157\\FW_8\\BTD\\btd_o_noisy%d_L%d',R(i),L(j));
save(temp,'Ubtd2');
plotme(yn')


a=load('H:\\OpenfMRI-ds157\\ds157\\sub-05\\func\\food.txt');
b=load('H:\\OpenfMRI-ds157\\ds157\\sub-05\\func\\non_food.txt');
TR=1.6;
a(:,1)=a(:,1)/TR;
b(:,1)=b(:,1)/TR;

t_cpd=y1n(:,6);
t_btd=squeeze(mean(yn_vis(:,5,:),1));
t_tpica=y11n(:,2);

for i=1:size(a,1)
%     max_f_cpd(i)=max(t_cpd(round(a(i,1))+1:round(a(i,1)+24/TR)+1,1));
%     max_nf_cpd(i)=max(t_cpd(round(b(i,1))+1:round(b(i,1)+24/TR)+1,1));
    max_f_btd(i)=max(t_btd(round(a(i,1))+1:round(a(i,1)+24/TR)+1,1));
    max_nf_btd(i)=max(t_btd(round(b(i,1))+1:round(b(i,1)+24/TR)+1,1));
%     max_f_tpica(i)=max(t_tpica(round(a(i,1))+1:round(a(i,1)+24/TR)+1,1));
%     max_nf_tpica(i)=max(t_tpica(round(b(i,1))+1:round(b(i,1)+24/TR)+1,1));
end
diff_cpd=mean(max_f_cpd)-mean(max_nf_cpd);
diff_btd=mean(max_f_btd)-mean(max_nf_btd);
diff_tpica=mean(max_f_tpica)-mean(max_nf_tpica);

tfood1=-4.9*ones(370,1);
tnfood1=4.9*ones(370,1);
tfood1([1:23,48:71,91:113,136:160,181:210,233:258,280:301,322:345],1)=-Uals{1,2}([1:23,48:71,91:113,136:160,181:210,233:258,280:301,322:345],4);
tnfood1([23:48,71:91,113:136,160:181,210:233,258:280,301:322,345:370],1)=Uals{1,2}([23:48,71:91,113:136,160:181,210:233,258:280,301:322,345:370],4);
%Semi blind
% Model the CPD of T, where the last columns of C are known
model = struct;
model.variables.a = permute(Utpica{1,1},[1,2,3,4,6,7,8,5]);
model.variables.b = permute(Utpica{1,2},[1,2,3,4,6,7,8,5]);
model.variables.b=model.variables.b(:,1:6);
model.variables.c = permute(Utpica{1,3},[1,2,3,4,6,7,8,5]);
model.factors.A = 'a';
model.factors.B = {'b',-tnfood1,tfood1}; % The third factor is partially known.
model.factors.C = 'c';
model.factorizations.myfac.data = dat1;
model.factorizations.myfac.cpd = {'A','B','C'};
sdf_check(model, 'print');
[sol,output]=sdf_nls(model,'Display',100);

Udriven{1,1}=sol.variables.a;
Udriven{1,2}=sol.factors.B;
Udriven{1,3}=sol.variables.c;


%PFAC2
for i=1:18
Pinit{1,i}=Uals1{1,2};
end
[A,H,C,P,fit]=parafac2(Sica,R,0,0,Uals1{1,1},eye(R),Uals1{1,3},Pinit);
clear a
for i=1:size(P,2)
    a(i,:,:)=P{1,i}*H;
    Upfac2{1,1}=A;
    Upfac2{1,2}=-squeeze(a(i,:,:));
    Upfac2{1,3}=C;
    Upfac2=cellfun(@(u,v)v*u,Upfac2,Vica,'UniformOutput',false); %%uncompressing results
    yy(i,:,:)=yn2';
    yn2=Upfac2{1,2};
    yn2=icatb_convertToZScores(yn2');
    yn2=icatb_filt_data(yn2',1,0.15);
    yn_vis(i,:,:)=yn2';
end



figure;plot(mean(yn_vis))
figure;plot(mean(yn_vis([1,2,3,4,5,7,8,10,11,12,13,15,16,18],:)));
timc=(squeeze(datall(:,:,14)))'*pinv(Upfac2{1,1})';
yn2=icatb_convertToZScores(timc(:,5)');
yn2=icatb_filt_data(yn2',1,0.15);

timc=(squeeze(datall(:,:,14)))'*pinv(Upfac2{1,1}(:,5))';
yn22=icatb_convertToZScores(timc');
yn22=icatb_filt_data(yn22',1,0.15);

figure; plot(-yn22); hold on; plot(-y1n(:,4),'r') 

Upfac2{1,1}(:,6)=(squeeze(datall(:,:,14)))*Upfac2{1,2}(:,5);