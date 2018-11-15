load('cpdtes_ica_parameter_info.mat')
load('H:\OpenfMRI-ds157\FW_8\CPD\components\cpdtes_ica_br1.mat')

Mask    = sesInfo.userInput.mask_ind;

map=compSet.ic;
coeff=zeros(size(map(1,:),2),1);

comp=reshape(map(6,:),d(1).hdr.dime.dim(2),d(1).hdr.dime.dim(3),[]);
[c,i]=max(comp(:));
[I1,I2,I3] = ind2sub(size(comp),i);

synth=(c-4)+8*rand([11*11*11,1]);
d3synt=zeros(size(comp));
d3synt(20:30,13:23,6:16)=reshape(synth,11,11,11);

d3synt(37:47,20:30,6:16)=reshape(synth,11,11,11);

d3synt=d3synt(:);
coeff(Mask)=d3synt(Mask);
map(1,:)=coeff;

a=load('H:\OpenfMRI-ds157\ds157\sub-05\func\food.txt');
b=load('H:\OpenfMRI-ds157\ds157\sub-05\func\non_food.txt');
TR=1.6;
a(:,1)=a(:,1)/TR+1;
b(:,1)=b(:,1)/TR+1;

tmap=zeros(size(compSet.it(1,:)));

tmap(1,[round(a(2,1)+3):round(a(2,1)+42/TR),round(a(4,1)+2):round(a(4,1)+42/TR),round(a(6,1)+3):round(a(6,1)+42/TR)])=3*ones(1,3*size(round(a(2,1)+3):round(a(2,1)+42/TR),2));
[hrf,p] = spm_hrf(1.6);
timec=conv(tmap,hrf);

comp1=map(1,:)'*timec(1,1:370);


% For non trilenar create shifted version 
timec2=circshift(conv(tmap,1.2*hrf),[0,7]);
comp2=map(1,:)'*timec2(1,1:370);

timec3=circshift(conv(tmap,0.7*hrf),[0,-7]);
comp3=map(1,:)'*timec3(1,1:370);

Save_Data2GUI(map,compSet.it',sesInfo);

%For semiblind
tmap=zeros(size(compSet.it(1,:)));
for i=1:size(a,1)
    tmap(1,round(a(i,1)):round(a(i,1)+24/TR))=ones(size([round(a(i,1)):round(a(i,1)+24/TR)]));
end
tfood=conv(tmap,hrf);
tmap=zeros(size(compSet.it(1,:)));
for i=1:size(a,1)
    tmap(1,round(b(i,1)):round(b(i,1)+24/TR))=ones(size([round(b(i,1)):round(b(i,1)+24/TR)]));
end
tnfood=conv(tmap,hrf);