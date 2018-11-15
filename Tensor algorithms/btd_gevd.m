function [ U_gevd1 ] = btd_gevd(d,R,L)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
U=cpd_gevd(reshape(d,[],size(d,3),size(d,4)),R);
 for k=1:R
   tempimg=reshape(U{1,1}(:,k),60,[]);
   [U1,S,V] = svd(tempimg);
   U_gevd1{1,k}{1,1}=U1(:,1:L)*S(1:L,1:L);
   U_gevd1{1,k}{1,2}=V(:,1:L);
   U_gevd1{1,k}{1,4}=U{1,3}(:,k);
   U_gevd1{1,k}{1,3}=U{1,2}(:,k);
   U_gevd1{1,k}{1,5}=diag(ones(1,L(1)));
 end

end

