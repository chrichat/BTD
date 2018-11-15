function [A B C err3]=els(X,A1,B1,C1,dA,dB,dC)
%
% computes parallel factor analysis
%
Ia=size(A1,1);
[Ib,r]=size(B1);
[Ic,r]=size(C1);
err3=chyba(X,A1,B1,C1);
T1=reshape(X,Ia,Ib*Ic)-A1*krb(C1,B1)';
T2=dA*krb(C1,B1)'+A1*(krb(dC,B1)+krb(C1,dB))';
T3=dA*(krb(dC,B1)+krb(C1,dB))'+A1*krb(dC,dB)';
T4=dA*krb(dC,dB)';
a1=-2*T1(:)'*T2(:);
a2=sum(T2(:).^2)-2*T1(:)'*T3(:);
a3=2*(T2(:)'*T3(:)-T1(:)'*T4(:));
a4=sum(T3(:).^2)+2*T2(:)'*T4(:);
a5=2*T3(:)'*T4(:);
a6=sum(T4(:).^2);
muall=roots((6:-1:1).*[a6 a5 a4 a3 a2 a1]);
muall=real(muall(abs(imag(muall))<eps));
[valmin,imin]=min(a1*muall+a2*muall.^2+a3*muall.^3+a4*muall.^4+a5*muall.^5+a6*muall.^6);
muvys=muall(imin); err3=valmin+err3;
A=A1+muvys*dA; B=B1+muvys*dB; C=C1+muvys*dC;
end

function ch=chyba(T,A,B,C)
%
n=size(T,1);
ch=sum((T(:)-reshape(A*krb(C,B)',n^3,1)).^2);
end
