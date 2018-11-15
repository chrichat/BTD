function X = psqrt(A,tol)

   % Produces A^(-.5) even if rank-problems

   [U,S,V] = svd(A,0);
   if min(size(S)) == 1
     S = S(1);
   else
     S = diag(S);
   end
   if (nargin == 1)
     tol = max(size(A)) * S(1) * eps;
   end
   r = sum(S > tol);
   if (r == 0)
     X = zeros(size(A'));
   else
     S = diag(ones(r,1)./sqrt(S(1:r)));
     X = V(:,1:r)*S*U(:,1:r)';
  end