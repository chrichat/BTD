function [Xwh, mu, invMat, whMat] = whiten(X,epsilon)
%function [X,mu,invMat] = whiten(X,epsilon)
if ~exist('epsilon','var')
    epsilon = 0.0001;
end

mu = mean(X); 
X = bsxfun(@minus, X, mu);
A = X'*X;
[V,D,notused] = svd(A);
whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
Xwh = X*whMat;  
invMat = pinv(whMat);

end
