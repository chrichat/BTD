load('fMRIsim-data_6060100_set1.mat')

C=5*rand(5,8);
A=reshape(SC,[],size(SC,3));
Ub{1,1}=A;
Ub{1,2}=W;

R=8;
L=20;
cols=ceil(sqrt(2*R));
rows=ceil(sqrt(2*R));

figure;
for i=1:size(SC,3)
  subplot(cols,rows,2*i-1)
  imagesc(reshape(A(:,i),60,60))
  subplot(cols,rows,2*i)
  plot(W(:,i))
end

X=A*kr(C,W)';
Es=randn(size(X));
c=norm(X,'fro')/(1.2*(norm(Es,'fro')));
X=X+c*Es;
X1=remmean(X);
data=reshape(X1,60*60,100,5);

% data1=data(:,:,1);
% data2=data(:,:,2);
% data3=data(:,:,3);
% data4=data(:,:,4);
% data5=data(:,:,5);
% [E1, D1]=pcamat(data1, 1, 8);
% [E2, D2]=pcamat(data2, 1, 8);
% [E3, D3]=pcamat(data3, 1, 8);
% [E4, D4]=pcamat(data4, 1, 8);
% [E5, D5]=pcamat(data5, 1, 8);
% 
% datas(:,:,1)=E1';
% datas(:,:,2)=E2';
% datas(:,:,3)=E3';
% datas(:,:,4)=E4';
% datas(:,:,5)=E5';
% [Wiva] = iva_laplace(datas);
% for i=1:size(C,1)
%     temp(:,:,i)=(Wiva(:,:,i))*datas(:,:,i);
%     tim(:,:,i)=data(:,:,i)'*pinv(temp(:,:,i));
% end
% k=3;
% figure
% for i=1:size(temp,1)
%   subplot(cols,rows,2*i-1)
%   imagesc(reshape(temp(i,:,k),60,60))
%   subplot(cols,rows,2*i)
%   plot(tim(:,i,k))
% end
% 
% temp1=mean(temp,3);
% tim1=mean(tim,3);
% figure;
% for i=1:size(temp,1)
%   subplot(cols,rows,2*i-1)
%   imagesc(reshape(temp1(i,:),60,60))
%   subplot(cols,rows,2*i)
%   plot(tim1(:,i))
% end
% Uiva{1,1}=temp1';
% Uiva{1,2}=tim1;
% [ corr_iva,scorr_iva,tcorr_iva ] = mean_corr( Uiva,Ub );
%  means_iva=mean(scorr_iva);
% meant_iva=mean(tcorr_iva);
%  
% if mod(cols,2)==1
%     cols=cols-1;
%     rows=rows+1;
% end

%Sparce ICA
lambda=10^-4;
epsilon=100;
tic
[icasigS] = ICA_EBM_Sparse(X,lambda,epsilon);
toc
% [icasig,A,B,ll,info]=icaML(data(:,:,1),R);
% [icasig,Rnew, Wnn] = fastica(X1','numOfIC',R); 
BS=data'*pinv(icasigS);
UtpicaS{1,1}=icasigS';
UtpicaS{1,2}=BS;
figure;
for i=1:size(icasigS,1)
  subplot(cols,rows,2*i-1)
  imagesc(reshape(icasigS(i,:),60,60))
  subplot(cols,rows,2*i)
  plot(-BS(:,i))
end
[ corr_gicaS,scorr_gicaS,tcorr_gicaS ] = mean_corr( UtpicaS,UbS );
means_gicaS=mean(scorr_gicaS);
meant_gicaS=mean(tcorr_gicaS);

% datass=permute(reshape(permute(datas,[2,1,3]),size(datas,2),[]),[2,1,3]);
% [Ess, Dss]=pcamat(datass', 1, 8);
[icasig,B,C1] = tpica(data,R,0);
% [icasig,A,B,ll,info]=icaML(data(:,:,1),R);
% [icasig,Rnew, Wnn] = fastica(X1','numOfIC',R); 
Utpica{1,1}=icasig';
Utpica{1,2}=B;
for i=1:size(C,1)
    Bnew(:,:,i)=data(:,:,i)'*pinv(icasig);
end
Bnewm=mean(Bnew,3)';
figure;
for i=1:size(icasig,1)
  subplot(cols,rows,2*i-1)
  imagesc(reshape(icasig(i,:),60,60))
  subplot(cols,rows,2*i)
  plot(-B(:,i))
end
[ corr_gica,scorr_gica,tcorr_gica ] = mean_corr( Utpica,Ub );
 means_gica=mean(scorr_gica);
meant_gica=mean(tcorr_gica);
 

 clear U1
 clear U2
 options.Compression = false;
 options.Algorithm=@cpd_nls;
 options.ExploitStructure=false;
 options.Refinement=false;
 options.Initialization = @cpd_gevd;
 options.Display=1;
 options.MaxIter=2000;
 [U1,output] = cpd(data,R,[196*3,196,3],options);
figure;
 for i=1:R
    subplot(cols,rows,2*i-1)
    imagesc(reshape(U1{1}(:,i),60,60))
    subplot(cols,rows,2*i)
    plot(U1{1,2}(:,i))
 end
 Ucpd=U1(1,1:2);
[ corr_cpd,scorr_cpd,tcorr_cpd ] = mean_corr( Ucpd,Ub );
 means_cpd=mean(scorr_cpd);
meant_cpd=mean(tcorr_cpd);
 

L=50;
N = [2 2]; % Order of the tensors
T4d=reshape(X1,60,60,[],size(data,3));
T=T4d;
U_mult=multiinit_lp(T,R,L,N);
Um2 = lp_nls(T,U_mult,[2 2]);
figure
for k=1:R
    tempimg=Um2{1,k}{1,1}*Um2{1,k}{1,2}';
    subplot(cols,rows,2*k-1)
    imagesc(tempimg)
    Ubtd{1,1}(:,k)=reshape(tempimg,1,[]);
    subplot(cols,rows,2*k)
    plot(Um2{1,k}{1,3})
    Ubtd{1,2}(:,k)=Um2{1,k}{1,3};
end
[ corr_btd,scorr_btd,tcorr_btd ] = mean_corr( Ubtd,Ub );
means_btd=mean(scorr_btd);
meant_btd=mean(tcorr_btd);