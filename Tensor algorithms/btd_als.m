function [ X,Y,B,C] = btd_als( T4d,X,Y,B,C,R,L,iter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
E=blkdiag(ones(1,L),ones(1,L),ones(1,L));
for j=1:iter
        W=((E'*(C'*C)*E).*(E'*(B'*B)*E).*(Y'*Y));
        X=tens2mat(T4d,1)*kr((C*E),(B*E),Y)*pinv(W);
        W=((E'*(C'*C)*E).*(E'*(B'*B)*E).*(X'*X));
        Y=tens2mat(T4d,2)*kr((C*E),(B*E),X)*pinv(W);
        W=E*((E'*(C'*C)*E).*(Y'*Y).*(X'*X))*E';
        B=tens2mat(T4d,3)*kr((C*E),Y,X)*E'*pinv(W);
        W=E*((E'*(B'*B)*E).*(Y'*Y).*(X'*X))*E';
        C=tens2mat(T4d,4)*kr((B*E),Y,X)*E'*pinv(W);
end
    for i=1:R
      A(:,i)=reshape(X(:,(i-1)*L+1:L*i)*Y(:,(i-1)*L+1:L*i)',[],1);
    end
 end

