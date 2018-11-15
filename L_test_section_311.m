spat=csvread('C:\Users\hatzi\Dropbox\matlab tools\Data\Helwig''s Data\Amat.csv');
B=csvread('C:\Users\hatzi\Dropbox\matlab tools\Data\Helwig''s Data\Bmat.csv');
brain_cond=csvread('C:\Users\hatzi\Dropbox\matlab tools\Data\Helwig''s Data\xybrain.csv');
SNR=[0.5, 0.15, 0.08];
patd=spat(:,3:6);
p=reshape(patd(:,4),60,50);
spatd=spat(:,1:2);

C=[3,4,5;2,3,4;2,2,3];

%Higher order source
a=reshape(patd(:,4),60,50);
a(52:55,12:15)=ones(4,4);
a(55:57,16:20)=ones(3,5);
a(41:45,31:35)=ones(5,5);
for i=1:3
    for k=1:30
        A=[spatd a(:)];
        Ub{1,1}=A;
        Ub{1,2}=B;
        aind=find(A);
        A(aind)=0.8+0.4.*rand(size(aind)); %set uniform activation from 0.8 to 1.2
        X=A*kr(C,B)';
        Es=randn(size(X));
        c2=norm(X,'fro')/(norm(SNR(i)*Es,'fro'));
        X=X+c2*Es;
        if i==3
            R=5;
        else
            R=3;
        end
        X1=remmean(X);
        data=reshape(X1,3000,180,3);
        options.Compression = false;
        options.Algorithm=@cpd_nls;
        options.ExploitStructure=false;
        options.Refinement=false;
        options.Initialization = @cpd_gevd;
        options.Display=10;
        options.MaxIter=2000;
        data=reshape(X1,3000,180,3);
        T4d=reshape(X1,60,50,[],size(data,3));
        T=T4d;
        for L=1:10
            U_mult=btd_gevd(T,R,L);
            Um2 = lp_nls(T,U_mult,[2 2]);
            for j=1:R
                tempimg(:,:,j)=Um2{1,j}{1,1}*Um2{1,j}{1,2}';
                Ubtd{1,1}(:,j)=reshape(tempimg(:,:,j),1,[]);
                Ubtd{1,1}(abs(Ubtd{1,1}(:,j))<0.25*max(abs(Ubtd{1,1}(:,j))),j)=0;
                Ubtd{1,2}(:,j)=Um2{1,j}{1,3};
            end
            [mcorr_btd(i,L),scorr_btd(:,L),tcorr_btd(:,L),ind(:,L),corrmatr_s_btd(:,:,L,k,i),corrmatr_t_btd(:,:,L,i)] = mean_corr(Ubtd,Ub);
            if L==3
                for j=1:size(Ub{1,1},2)
                    if abs(min(Ubtd{1,2}(:,ind(j,L))))>abs(max(Ubtd{1,2}(:,ind(j,L))))
                        Ubtd{1,2}(:,ind(j,L))=-(Ubtd{1,2}(:,ind(j,L)));
                    end
                    if abs(min(min(tempimg(:,:,ind(j,L)))))>abs(max(max(tempimg(:,:,ind(j,L)))))
                        tempimg(:,:,ind(j,L))=-tempimg(:,:,ind(j,L));
                    end
                end
            end
        end
        for j=1:L
            mean_L(j,k,i)=mean(corrmatr_s_btd(3,3,j,k,:));
        end
    end
end