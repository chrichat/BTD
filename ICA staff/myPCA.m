function res = myPCA(X)
[D, N] = size(X);
m = mean(X, 2);
X = X - repmat(m, 1, N);
[res, s, v] = svd(X,'econ');
% newu=e(:,1:k);
% res=newu'*X;
end