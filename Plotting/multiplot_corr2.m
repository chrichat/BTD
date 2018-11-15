function multiplot_corr2( cpd2s, cpd2t, btd2s, btd2t, autocors,autocort )

%%% Input- 
% cpds2=space autocorrelation of cpd 2
% cpd2t=time autocorrelation of cpd2 
% btd2s=space autocorrelation of btd2 
% btd2t=time autocorrelation of btd2
%autocors=autocorrelation space
%autocort=autocorrelation time

cpd2s=permute(cpd2s,[3,2,1]);
cpd2t=permute(cpd2t,[3,2,1]);
btd2s=permute(btd2s,[3,2,1]);
btd2t=permute(btd2t,[3,2,1]);

autcorr_cpd2=zeros(size(cpd2s,1),size(cpd2s,2));
autcorr_btd2=zeros(size(cpd2s,1),size(cpd2s,2));
for k=1:size(cpd2s,3)
    for i=1:size(cpd2s,2)
        for j=1:size(cpd2s,1)
            if i==j
                corr_cpd2(i,k)=(abs(cpd2s(i,j,k)-autocors(i,j))+abs(cpd2t(i,j,k)-autocort(i,j)))/2;
                corr_btd2(i,k)=(abs(btd2s(i,j,k)-autocors(i,j))+abs(btd2t(i,j,k)-autocort(i,j)))/2;
            else
                autcorr_cpd2(i,k)=autcorr_cpd2(i,k)+(abs(cpd2s(i,j,k)-autocors(i,j))/2+abs(cpd2t(i,j,k)-autocort(i,j))/2)/2;
                autcorr_btd2(i,k)=autcorr_btd2(i,k)+(abs(btd2s(i,j,k)-autocors(i,j))/2+abs(btd2t(i,j,k)-autocort(i,j))/2)/2;
            end
            
        end
    end
end

figure
xx=[1 , 1.1, 1.2, 1.4 ,1.5, 1.6, 1.8, 1.9, 2.0];
scatter(xx,corr_cpd2(:),200,[1 0 0 ],'filled')
hold on
scatter(xx,corr_btd2(:),200,[0.4 0.8 0.4],'filled')
scatter(xx,autcorr_cpd2(:),200,[1 0 0 ],'d','filled')
scatter(xx,autcorr_btd2(:),200,[0.4 0.8 0.4],'d','filled')
h=legend({'PARAFAC2;Principal CCD', 'BTD2;Principal CCD','PARAFAC2;Mean crosstalk CCD','BTD2;Mean crosstalk CCD'});
set(h,'FontSize',16)
legendmarkeradjust(20)
set(gca,'Xtick',[0.9, 1 , 1.1, 1.2, 1.3, 1.4 ,1.5, 1.6, 1.7 ,1.8, 1.9, 2.0, 2.1])
set(gca,'Xticklabel',{'','Source 1','Source 2','Source 3','','Source 1','Source 2','Source 3','','Source 1','Source 2','Source 3',''})
set(gca,'FontSize',16,'FontWeight','Bold')


figure
xxline=[1,2,3,5,6,7,9,10,11];
bar(xxline-0.2,1-(corr_btd2(:)),0.10,'FaceColor',[0.345, 0.549, 0.450]);
hold on
scatter(xxline-0.2,1-autcorr_btd2(:),150,[0.345, 0.549, 0.450],'d','filled','MarkerEdgeColor','k');
bar(xxline,1-(corr_cpd2(:)),0.10,'FaceColor',[0.949, 0.890, 0.580]);
scatter(xxline,1-autcorr_cpd2(:),150,[0.949, 0.890, 0.580],'d','filled','MarkerEdgeColor','k');
h=legend({'BTD2;Principal','BTD2;Crosstalk','CPD2;Principal CCD','CPD2;Crosstalk'});
set(h,'FontSize',16)
legendmarkeradjust(20)
set(gca,'Xtick',[1:11])
set(gca,'Xticklabel',{'Source 1','Source 2','Source 3','','Source 1','Source 2','Source 3','','Source 1','Source 2','Source 3',''})
set(gca,'FontSize',16,'FontWeight','Bold')
end

