function [A,B,C,iter]=palsE(T,r,NN,A,B,C)
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
  A0=A; B0=B; C0=C;  
  ii=rand(1,r);
  i1=ii<1/3; i3=ii>2/3; i2=(i1+i3)==0;
  i1=find(i1); i2=find(i2); i3=find(i3);
  [A B C err1]=ajdr(T,A,B,C,i1,i2,i3);
  [A B C err]=ajdr(T,A,B,C,i2,i3,i1);
  [A B C err2]=ajdr(T,A,B,C,i3,i1,i2);
  [A B C err]=els(T,A0,B0,C0,A-A0,B-B0,C-C0);
  iter(i)=err;
  [i err1 err2 err]
end    
semilogy(iter,'k')
err
end