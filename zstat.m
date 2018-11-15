function L=zstat(Um2,az)
 % Code for Heuristic 1
 %
 %- Input 
 %---Um2 is the output of the lp_nls algorithm. It shall be a cell
 %type structure with R different cells, each one containing one source.
 %Each one of the different R cells must contain 5 cells; 2 cells representing the
 %two spatial domains with dimensions (I_x x R) and (I_yz x R)
 %respectively, one cell for the time domain (I_t x R), one for the subject
 %domain (I_s x R) and the last cell a diagonal matrix ( R x R). For more
 %details have a look at funcion lp_nls
 %---az. The p value used for testing the null hypothesis. Typical values
 %for P-value is 0.1, 0.05 and 0.01. We recommend the use of 0.05 in cases
 %with high instance of noise az=0.01 shall be used.
 
 
for j=1:size(Um2,2)
a=Um2{1,j}{1,1};
b=Um2{1,j}{1,2};
Lin=size(b,2);
for i=1:Lin
    c(:,:,i)=a(:,i)*b(:,i)';
end
ct(:,:,j)=reshape(c,[],Lin);
end

p=az/(size(ct,1)*size(ct,2)*size(ct,3)); %Bonferonni corrected p value
z=norminv(p);

for j=1:size(Um2,2)
cnew2=ct(:,:,j)-mean(ct(:));
cnew3=cnew2./std(ct(:));

bb=(abs(cnew3)-abs(z)+abs(abs(cnew3)-abs(z)));
bb(bb>0)=ones(size(find(bb>0)));
figure;
imagesc(bb)
[~, colIdcs] = find(bb~=0);
c=unique(colIdcs); 
Li(j)=size(c,1);
end
L=max(Li);

end

