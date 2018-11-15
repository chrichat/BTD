function [A,B,C,iter]=pals(T,r,NN,A,B,C)
%
[Ia,Ib,Ic]=size(T);
if nargin<3
    NN=200;
end  
if nargin<4
   A=randn(Ia,r);
   B=randn(Ib,r);
   C=randn(Ic,r);
end   
a1=sum(A.*conj(A)); b1=sum(B.*conj(B)); c1=sum(C.*conj(C));
p1=(a1.*b1.*c1).^(1/3);
A=A.*repmat(sqrt(p1./a1),Ia,1); %% normalize the initial matrices
B=B.*repmat(sqrt(p1./b1),Ib,1);
C=C.*repmat(sqrt(p1./c1),Ic,1);
iter=zeros(1,NN);
for i=1:NN
  ii=rand(1,r);
  i1=ii<1/3; i3=ii>2/3; i2=(i1+i3)==0;
  i1=find(i1); i2=find(i2); i3=find(i3);
  [A B C err]=ajdr(T,A,B,C,i1,i2,i3);
  [A B C err]=ajdr(T,A,B,C,i2,i3,i1);
  [A B C err]=ajdr(T,A,B,C,i3,i1,i2);
  iter(i)=err;
end    
semilogy(iter)
err
end
%
function [A,B,C,err]=ajdr(T,A,B,C,i1,i2,i3)
%
% computes Jacobian for symmetric CP decomposition of the nasob(n) tensor
%
[n,r]=size(A);
ic1=length(i1); ic2=length(i2); ic3=length(i3);
% J=zeros(n^3,n*r);
% if ic1>0
%    J(:,1:ic1*n)=kron(krb(C(:,i1),B(:,i1)),eye(n));
% end
% if ic2>0
%   for ir=1:ic2
%     ipoc=(ic1+ir-1)*n; 
%     J(:,ipoc+1:ipoc+n)=kron(C(:,i2(ir)),kron(eye(n),A(:,i2(ir))));
%   end
% end
% ic3=length(i3);
% if ic3>0
%   for ir=1:ic3
%     ipoc=(ic1+ic2+ir-1)*n; 
%     J(:,ipoc+1:ipoc+n)=kron(eye(n),kron(B(:,i3(ir)),A(:,i3(ir))));
%   end
% end
% H=J'*J;
I=eye(n);
HAA=kron((B(:,i1)'*B(:,i1)).*(C(:,i1)'*C(:,i1)),I);
HBB=kron((A(:,i2)'*A(:,i2)).*(C(:,i2)'*C(:,i2)),I);
HCC=kron((B(:,i3)'*B(:,i3)).*(A(:,i3)'*A(:,i3)),I);
%HAC=zeros(ic1*n,ic3*n); HBC=zeros(ic2*n,ic3*n); 
HBA=vec(A(:,i2))*vec(B(:,i1))';
aux1=reshape(HBA,n,ic2,n,ic1); aux2=permute(aux1,[3 2 1 4]);
if ic1*ic2>0
   HBA=reshape(aux2,n*ic2,n*ic1).*kron(C(:,i2)'*C(:,i1),ones(n,n));
else
   HBA=reshape(aux2,n*ic2,n*ic1);
end   
HCA=vec(A(:,i3))*vec(C(:,i1))';
aux1=reshape(HCA,n,ic3,n,ic1); aux2=permute(aux1,[3 2 1 4]);
if ic1*ic3>0
    HCA=reshape(aux2,n*ic3,n*ic1).*kron(B(:,i3)'*B(:,i1),ones(n,n));
else
    HCA=reshape(aux2,n*ic3,n*ic1);
end    
HCB=vec(B(:,i3))*vec(C(:,i2))';
aux1=reshape(HCB,n,ic3,n,ic2); aux2=permute(aux1,[3 2 1 4]);
if ic3*ic2>0
    HCB=reshape(aux2,n*ic3,n*ic2).*kron(A(:,i3)'*A(:,i2),ones(n,n));
else
    HCB=reshape(aux2,n*ic3,n*ic2);
end    
HH=[HAA HBA' HCA'; HBA HBB HCB'; HCA HCB HCC];
%g=J'*T(:);
g1=reshape(T,n,n^2)*krb(C(:,i1),B(:,i1));
g2=reshape(permute(T,[2,1,3]),n,n^2)*krb(C(:,i2),A(:,i2));
g3=reshape(T,n^2,n)'*krb(B(:,i3),A(:,i3));
g=[g1(:); g2(:); g3(:)];
%ff=H\g;
ff=HH\g;
%err0=sum((T(:)-reshape(A*krb(C,B)',n^3,1)).^2)
A(:,i1)=reshape(ff(1:ic1*n),n,ic1);
B(:,i2)=reshape(ff(ic1*n+1:(ic1+ic2)*n),n,ic2);
C(:,i3)=reshape(ff((ic1+ic2)*n+1:r*n),n,ic3);
err=sum((T(:)-reshape(A*krb(C,B)',n^3,1)).^2);
end

