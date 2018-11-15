function [ A,B,C,fit] = btd2(data,X,Y,H,C,R,L,maxit,iter )
%UNTITLED BTD2 using als
K=size(data,4);
d=reshape(data,[],size(data,3),size(data,4));
for i=1:R
    A(:,i)=reshape(X(:,(i-1)*L+1:L*i)*Y(:,(i-1)*L+1:L*i)',[],1);
end
for z=1:maxit
    clear T4
    clear P
    for k = 1:K
        Qk= d(:,:,k)'*(A*diag(C(k,:))*H');
        P{k}= Qk*psqrt(Qk'*Qk);
        T4(:,:,k) = d(:,:,k)*P{k};
    end
    d4=reshape(T4,size(data,1),size(data,2),size(T4,2),[]);
    [ X,Y,H,C] = btd_als(d4,X,Y,H,C,R,L,iter);
    for i=1:R
    A(:,i)=reshape(X(:,(i-1)*L+1:L*i)*Y(:,(i-1)*L+1:L*i)',[],1);
    end

    fit=0;
    for k = 1:K
        B=P{k}*H;
        M   = A*diag(C(k,:))*(P{k}*H)';
        fit = fit + sum(sum(abs (d(:,:,k) - M ).^2));
    end
    sprintf('The fit is %d', fit)
end
end


