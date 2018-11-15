function [ mcorr,scorr,tcorr,ind,corrmatr_s,corrmatr_t] = mean_corr( U1,U2 )
%Function which computes the mean correlation, the maximum correlation in
%space and time and provides reoriented the correlation martix.

tempcorr=abs(corr(U2{1,1},U1{1,1}));
[scorr,ind]=max(tempcorr,[],2);
[~,id_un] = unique(ind,'rows','stable'); % Find unique values
dup_id=setdiff(1:size(ind,1),id_un); %Find replicated sources and then replace them
% Check if some of the correlations are with the same map
while size(dup_id)~=0
    inds=find(ind==ind(dup_id(1)));
    [~,tempid]=min(scorr(inds));
    tempcorr(inds(tempid),ind(inds(tempid)))=0;
    [scorr,ind]=max(tempcorr,[],2);
    [~,id_un] = unique(ind,'rows','stable'); % Find unique values
    dup_id=setdiff(1:size(ind,1),id_un); %Find replicated sources and then replace them
end
corrmatr_s=[];
for j=1:size(ind,1)
  corrmatr_s=[corrmatr_s;tempcorr(:,ind(j))'];  
end
     
tempcorr=abs(corr(U2{1,2},U1{1,2}));
[tcorr,ind]=max(tempcorr,[],2);
[~,id_un] = unique(ind,'rows','stable'); % Find unique values
dup_id=setdiff(1:size(ind,1),id_un); %Find replicated sources and then replace them
% Check if some of the correlations are with the same map
while size(dup_id)~=0
    inds=find(ind==ind(dup_id(1)));
    [~,tempid]=min(tcorr(inds));
    tempcorr(inds(tempid),ind(inds(tempid)))=0;
    [tcorr,ind]=max(tempcorr,[],2);
    [~,id_un] = unique(ind,'rows','stable'); % Find unique values
    dup_id=setdiff(1:size(ind,1),id_un); %Find replicated sources and then replace them
end
mcorr=mean([mean(scorr),mean(tcorr)]);    
corrmatr_t=[];
for j=1:size(ind,1)
  corrmatr_t=[corrmatr_t;tempcorr(:,ind(j))'];  
end
end

