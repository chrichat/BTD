function multiplot_corr5( tpicas, tpicat, cpds, cpdt, btds, btdt, sicas,sicat, cpaicas,cpaicat, autocors,autocort )

%%% Input- 
% tpicas=space autocorrelation of tpica 
% tpicat=time autocorrelation of tpica 
% cpds=space autocorrelation of cpd 
% cpdt=time autocorrelation of cpd 
% btds=space autocorrelation of btd 
% btdt=time autocorrelation of btd
% cpaicas=space autocorrelation of cpa-ica 
% cpaicat=time autocorrelation of cpa-ica
% sica=space autocorrelation of sica 
% sicat=time autocorrelation of sica
%autocors=autocorrelation space
%autocort=autocorrelation time

autcorr_tpica=zeros(size(tpicas,1),size(tpicas,2));
autcorr_cpd=zeros(size(tpicas,1),size(tpicas,2));
autcorr_btd=zeros(size(tpicas,1),size(tpicas,2));
autcorr_sica=zeros(size(tpicas,1),size(tpicas,2));
autcorr_cpaica=zeros(size(tpicas,1),size(tpicas,2));

for k=1:size(tpicas,3)
    for i=1:size(tpicas,2)
        for j=1:size(tpicas,1)
            if i==j
                corr_tpica(i,k)=(abs(tpicas(i,j,k)-autocors(i,j))+abs(tpicat(i,j,k)-autocort(i,j)))/2;
                corr_cpd(i,k)=(abs(cpds(i,j,k)-autocors(i,j))+abs(cpdt(i,j,k)-autocort(i,j)))/2;
                corr_btd(i,k)=(abs(btds(i,j,k)-autocors(i,j))+abs(btdt(i,j,k)-autocort(i,j)))/2;
                corr_sica(i,k)=(abs(sicas(i,j,k)-autocors(i,j))+abs(sicat(i,j,k)-autocort(i,j)))/2;
                corr_cpaica(i,k)=(abs(cpaicas(i,j,k)-autocors(i,j))+abs(cpaicat(i,j,k)-autocort(i,j)))/2;

            else
                autcorr_tpica(i,k)=autcorr_tpica(i,k)+(abs(tpicas(i,j,k)-autocors(i,j))/2+abs(tpicat(i,j,k)-autocort(i,j))/2)/2;
                autcorr_cpd(i,k)=autcorr_cpd(i,k)+(abs(cpds(i,j,k)-autocors(i,j))/2+abs(cpdt(i,j,k)-autocort(i,j))/2)/2;
                autcorr_btd(i,k)=autcorr_btd(i,k)+(abs(btds(i,j,k)-autocors(i,j))/2+abs(btdt(i,j,k)-autocort(i,j))/2)/2;
                autcorr_sica(i,k)=autcorr_sica(i,k)+(abs(sicas(i,j,k)-autocors(i,j))/2+abs(sicat(i,j,k)-autocort(i,j))/2)/2;
                autcorr_cpaica(i,k)=autcorr_cpaica(i,k)+(abs(cpaicas(i,j,k)-autocors(i,j))/2+abs(cpaicat(i,j,k)-autocort(i,j))/2)/2;

            end
            
        end
    end
end


figure
xxline=[0.5,1.2,1.9];
bar(xxline-0.2,1-mean(autcorr_btd),0.08,'FaceColor',[0.345, 0.549, 0.450]);
hold on
errorbar(xxline-0.2,1-mean(autcorr_btd),mean(autcorr_btd)-max(autcorr_btd),min(autcorr_btd)-mean(autcorr_btd),'o','MarkerSize',.1,...
    'MarkerEdgeColor',[0.345, 0.549, 0.450],'Color','black','LineWidth', 1)
bar(xxline-0.1,1-mean(corr_cpd),0.08,'FaceColor',[0.949, 0.890, 0.580])
errorbar(xxline-0.1,1-mean(corr_cpd),mean(corr_cpd)-max(corr_cpd),min(corr_cpd)-mean(corr_cpd),'o','MarkerSize',.1,...
    'MarkerEdgeColor',[0.949, 0.890, 0.580],'Color','black','LineWidth', 1)
bar(xxline,1-mean(corr_sica),0.08,'FaceColor',[0.850, 0.392, 0.349])
errorbar(xxline,1-mean(corr_sica),mean(corr_sica)-max(corr_sica),min(corr_sica)-mean(corr_sica),'o','MarkerSize',.1,...
    'MarkerEdgeColor',[0.850, 0.392, 0.349],'Color','black','LineWidth', 1)
bar(xxline+0.1,1-mean(corr_cpaica),0.08,'FaceColor',[0.549, 0.274, 0.274])
errorbar(xxline+0.1,1-mean(corr_cpaica),mean(corr_cpaica)-max(corr_cpaica),min(corr_cpaica)-mean(corr_cpaica),'o','MarkerSize',.1,...
    'MarkerEdgeColor',[0.549, 0.274, 0.274],'Color','black','LineWidth', 1)
bar(xxline+0.2,1-mean(corr_tpica),0.08,'FaceColor',[0.949, 0.682, 0.447])
errorbar(xxline+0.2,1-mean(corr_tpica),mean(corr_tpica)-max(corr_tpica),min(corr_tpica)-mean(corr_tpica),'o','MarkerSize',.1,...
    'MarkerEdgeColor',[0.949, 0.682, 0.447],'Color','black','LineWidth', 1)

set(gca,'Xtick',[0.5,0.8,1.2,1.6,1.9])
set(gca,'Xticklabel',{'CNR=0.5','','CNR=0.15','','CNR=0.08'})

figure
xxline=[1, 2, 3, 5 ,6 ,7 ,9, 10 ,11];
scatter(xxline-0.3,corr_btd(:),200,[0.345, 0.549, 0.450],'filled');
hold on
scatter(xxline-0.3,autcorr_btd(:),200,[0.345, 0.549, 0.450],'d','filled');

scatter(xxline-0.15,corr_cpd(:),200,[0.949, 0.890, 0.580],'filled');
scatter(xxline-0.15,autcorr_cpd(:),200,[0.949, 0.890, 0.580],'d','filled');

scatter(xxline,corr_sica(:),200,[0.345, 0.549, 0.450],'filled');
scatter(xxline,autcorr_sica(:),200,[0.345, 0.549, 0.450],'d','filled');

scatter(xxline+0.15,corr_cpaica(:),200,[0.549, 0.274, 0.274],'filled');
scatter(xxline+0.15,autcorr_cpaica(:),200,[0.549, 0.274, 0.274],'d','filled');

scatter(xxline+0.3,corr_tpica(:),200,[0.949, 0.682, 0.447],'filled');
scatter(xxline+0.3,autcorr_tpica(:),200,[0.949, 0.682, 0.447],'d','filled');


figure
xxline=[1:9];
scatter(xxline-0.3,1-corr_btd(:),200,[0.345, 0.549, 0.450],'filled');
hold on
scatter(xxline-0.3,1-autcorr_btd(:),200,[0.345, 0.549, 0.450],'d','filled');

scatter(xxline-0.15,1-corr_cpd(:),200,[0.949, 0.890, 0.580],'filled');
scatter(xxline-0.15,1-autcorr_cpd(:),200,[0.949, 0.890, 0.580],'d','filled');

scatter(xxline,1-corr_sica(:),200,[0.345, 0.549, 0.450],'filled');
scatter(xxline,1-autcorr_sica(:),200,[0.345, 0.549, 0.450],'d','filled');

scatter(xxline+0.15,1-corr_cpaica(:),200,[0.549, 0.274, 0.274],'filled');
scatter(xxline+0.15,1-autcorr_cpaica(:),200,[0.549, 0.274, 0.274],'d','filled');

scatter(xxline+0.3,1-corr_tpica(:),200,[0.949, 0.682, 0.447],'filled');
scatter(xxline+0.3,1-autcorr_tpica(:),200,[0.949, 0.682, 0.447],'d','filled');

figure
xxline=[1, 2.5, 4, 6, 7.5, 9, 11, 12.5, 14];
bar(xxline-0.5,1-(corr_btd(:)),0.10,'FaceColor',[0.345, 0.549, 0.450]);
hold on
% h=bar(xxline-0.41,1-(autcorr_btd(:)),0.06,'FaceColor',[0.345, 0.549, 0.450]);
% hatchfill2(h,'single','HatchAngle',45);
scatter(xxline-0.5,1-autcorr_btd(:),150,[0.345, 0.549, 0.450],'d','filled','MarkerEdgeColor','k');

bar(xxline-0.25,1-(corr_cpd(:)),0.10,'FaceColor',[0.949, 0.890, 0.580]);
% h=bar(xxline-0.15,1-(autcorr_cpd(:)),0.06,'FaceColor',[0.949, 0.890, 0.580]);
% hatchfill2(h,'single','HatchAngle',45);
scatter(xxline-0.25,1-autcorr_cpd(:),150,[0.949, 0.890, 0.580],'d','filled','MarkerEdgeColor','k');

bar(xxline,1-(corr_sica(:)),0.10,'FaceColor',[0.850, 0.392, 0.349]);
% h=bar(xxline+0.16,1-(autcorr_sica(:)),0.06,'FaceColor',[0.850, 0.392, 0.349]);
% hatchfill2(h,'single','HatchAngle',45);
scatter(xxline,1-autcorr_sica(:),150,[0.850, 0.392, 0.349],'d','filled','MarkerEdgeColor','k');

bar(xxline+0.25,1-(corr_cpaica(:)),0.10,'FaceColor',[0.549, 0.274, 0.274]);
% h=bar(xxline+0.41,1-(autcorr_cpaica(:)),0.06,'FaceColor',[0.549, 0.274, 0.274]);
% hatchfill2(h,'single','HatchAngle',45);
scatter(xxline+0.25,1-autcorr_cpaica(:),150,[0.549, 0.274, 0.274],'d','filled','MarkerEdgeColor','k');

bar(xxline+0.5,1-(corr_tpica(:)),0.10,'FaceColor',[0.949, 0.682, 0.447]);
% h=bar(xxline+0.66,1-(autcorr_tpica(:)),0.06,'FaceColor',[0.949, 0.682, 0.447]);
% hatchfill2(h,'single','HatchAngle',45);
scatter(xxline+0.5,1-autcorr_tpica(:),150,[0.949, 0.682, 0.447],'d','filled','MarkerEdgeColor','k');


% figure
% herrorbar(xxline,1-mean(corr_cpd),std(corr_btd)/2)
% hold on
% herrorbar(xxline,1-mean(corr_btd),std(corr_btd)/2)
% herrorbar(xxline,1-mean(corr_tpica),std(corr_tpica)/2)
% herrorbar(xxline,1-mean(corr_sica),std(corr_sica)/2)
% herrorbar(xxline,1-mean(corr_cpaica),std(corr_cpaica)/2)
% set(gca,'Xtick',[0.08,0.15,0.5])
% set(gca,'Xticklabel',{'0.08','0.15','0.5'})
% set(gca, 'XDir','reverse')

end

