% spat=csvread('C:\Users\hatzi\Dropbox\matlab tools\Data\Helwig''s Data\Amat.csv');
% B=csvread('C:\Users\hatzi\Dropbox\matlab tools\Data\Helwig''s Data\Bmat.csv');
% brain_cond=csvread('C:\Users\hatzi\Dropbox\matlab tools\Data\Helwig''s Data\xybrain.csv');

SNR=[0.5, 0.15, 0.08]; %Different levels of noise to be tested
patd=spat(:,3:6);
p=reshape(patd(:,4),60,50);
spatd=spat(:,1:2); %use only the patd with the highest overlap

C=[3,4,5;2,3,4;2,2,3]; %Set values of the subject amplitudes for each source

%Higher order source - Change source 3 to a source with higher order
a=reshape(patd(:,4),60,50);
a(52:55,12:15)=ones(4,4);
a(55:57,16:20)=ones(3,5);
a(41:45,31:35)=ones(5,5);

%Mix the sources
A=[spatd a(:)];
Ub{1,1}=A;
Ub{1,2}=B;
aind=find(A);
A(aind)=0.8+0.4.*rand(size(aind)); %set uniform activation from 0.8 to 1.2
X=A*kr(C,B)';
    
figure('rend','painters','pos',[10 10 1000 700])
for i=1:size(SNR,2)
    
    %Add noise
    Es=randn(size(X));
    c2=norm(X,'fro')/(norm(SNR(i)*Es,'fro'));
    X=X+c2*Es;
        if i==3
            R=5;
        else
            R=3;
        end
     X1=remmean(X');
     
%      % Test sparse EBM
%      tic
%      [E, D]=pcamat(X1, 1,3);
%      [whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
% 						     (X1, E, D);
%      lambda=10^-2;
%      epsilon=10^-1;    
%      [W] = ICA_EBM_Sparse(whitesig,lambda,epsilon);
%      toc
%      UtpicaS{1,1} = W * whitesig ;
%      plotme(UtpicaS{1,1},60,50)
%      %
%      
     
    %TPICA with INFOMAX
     tic
     bb=icaML(X1,3);
     Utpica{1,1}=bb';
     temp=(X'*pinv(bb));
     toc
     Utpica{1,2}=temp(1:180,:);
     [mcorr,scorr,tcorr,ind,corrmatr_s_ica(:,:,i),corrmatr_t_ica(:,:,i)] = mean_corr(Utpica,Ub);
     for j=1:size(Ub{1,1},2)
         if abs(min(Utpica{1,2}(:,ind(j))))>abs(max(Utpica{1,2}(:,ind(j))))
            Utpica{1,2}(:,ind(j))=-(Utpica{1,2}(:,ind(j)));
         end
         if abs(min(Utpica{1,1}(:,ind(j))))>abs(max(Utpica{1,1}(:,ind(j))))
            Utpica{1,1}(:,ind(j))=-Utpica{1,1}(:,ind(j));
         end
         subplot(size(SNR,2)*3,6,(i-1)*18+(j-1)*6+1)
         imagesc(flipud(reshape(Utpica{1,1}(:,ind(j)),60,50)))
         colorbar
         hold on
         plot(brain_cond(:,1),brain_cond(:,2),'w','LineWidth',3)
         set(gca,'YDir','normal')
         hold off
         subplot(size(SNR,2)*3,6,(i-1)*18+(j-1)*6+2)
         plot(Utpica{1,2}(:,ind(j)))
         axis tight
     end
     clear U1
     clear U2
     
     %CPD
     X1=remmean(X);
     data=reshape(X1,3000,180,3);
     options.Compression = false;
     options.Algorithm=@cpd_nls;
     options.ExploitStructure=false;
     options.Refinement=false;
     options.Initialization = @cpd_gevd;
     options.Display=100;
     options.MaxIter=2000;
     tic
     [U1,output] = cpd(data,R,options);
     toc
     Ucpd=U1(1,1:2);
     [mcorr_cpd,scorr_cpd,tcorr_cpd,ind,corrmatr_s_cpd(:,:,i),corrmatr_t_cpd(:,:,i)] = mean_corr(Ucpd,Ub);
     for j=1:size(Ub{1,1},2)
         if abs(min(Ucpd{1,2}(:,ind(j))))>abs(max(Ucpd{1,2}(:,ind(j))))
            Ucpd{1,2}(:,ind(j))=-(Ucpd{1,2}(:,ind(j)));
         end
         if abs(min(Ucpd{1,1}(:,ind(j))))>abs(max(Ucpd{1,1}(:,ind(j))))
            Ucpd{1,1}(:,ind(j))=-Ucpd{1,1}(:,ind(j));
         end
         subplot(size(SNR,2)*3,6,(i-1)*18+(j-1)*6+3)
         imagesc(flipud(reshape(Ucpd{1,1}(:,ind(j)),60,50)))
         colorbar
         hold on
         plot(brain_cond(:,1),brain_cond(:,2),'w','LineWidth',3)
         set(gca,'YDir','normal')
         hold off
         subplot(size(SNR,2)*3,6,(i-1)*18+(j-1)*6+4)
         plot(Ucpd{1,2}(:,ind(j)))
         axis tight
    end
    data=reshape(X1,3000,180,3);
    
    %BTD
    T4d=reshape(X1,60,50,[],size(data,3));
    T=T4d;
    for L=1:10
        tic
        U_mult=btd_gevd(T,R,L);
        Um2 = lp_nls(T,U_mult,[2 2]);
        toc
        for j=1:R    
            tempimg(:,:,j)=Um2{1,j}{1,1}*Um2{1,j}{1,2}';
            Ubtd{1,1}(:,j)=reshape(tempimg(:,:,j),1,[]);
            Ubtd{1,1}(abs(Ubtd{1,1}(:,j))<0.25*max(abs(Ubtd{1,1}(:,j))),j)=0;
            Ubtd{1,2}(:,j)=Um2{1,j}{1,3};
        end
        [mcorr_btd(L),scorr_btd(:,L),tcorr_btd(:,L),ind(:,L),corrmatr_s_btd(:,:,L,i),corrmatr_t_btd(:,:,L,i)] = mean_corr(Ubtd,Ub);
        if L==3
            for j=1:size(Ub{1,1},2)
                if abs(min(Ubtd{1,2}(:,ind(j,L))))>abs(max(Ubtd{1,2}(:,ind(j,L))))
                    Ubtd{1,2}(:,ind(j,L))=-(Ubtd{1,2}(:,ind(j,L)));
                end
                if abs(min(min(tempimg(:,:,ind(j,L)))))>abs(max(max(tempimg(:,:,ind(j,L)))))
                    tempimg(:,:,ind(j,L))=-tempimg(:,:,ind(j,L));
                end
                subplot(size(SNR,2)*3,6,(i-1)*18+(j-1)*6+5)
                imagesc(flipud(tempimg(:,:,ind(j,L))))
                colorbar
                hold on
                plot(brain_cond(:,1),brain_cond(:,2),'w','LineWidth',3)
                set(gca,'YDir','normal')
                hold off
                subplot(size(SNR,2)*3,6,(i-1)*18+(j-1)*6+6)
                plot(Ubtd{1,2}(:,ind(j,L)))
                axis tight
            end
        end
    end
end
 