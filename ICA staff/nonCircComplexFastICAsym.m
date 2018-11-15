%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference
% Mike Novey and T. Adali, "On Extending the complex FastICA algorithm
% to noncircular sources" IEEE Trans. Signal Processing, vol. 56, no. 5, pp. 2148-2154, May 2008.
%
% Performs symmetric orthogonalization of the
% non-circular complex FastICA algorithm where
% xold is the mixtures and typeStr = 'log', 'kurt', or 'sqrt' is
% the nonlinearity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Ahat, shat] = nonCircComplexFastICAsym(xold,typeStr)
 type = 0;
 if strcmp(typeStr,'log') == 1
     type = 1;
 elseif strcmp(typeStr,'kurt') == 1
     type = 2;
 elseif strcmp(typeStr,'sqrt') == 1
     type = 3;
 end
 
 tol = 1e-5;
 a2=.05;
  defl = 1; % components are estimated one by one in a deflationary manner; set this to 0 if you want them all estimated simultaneously
 maxcounter = 50;
  [n,m] = size(xold);
% Whitening of s:
yyy = zeros(1,m);
[Ex, Dx] = eig(cov(xold'));
Q = mtimes(sqrt(inv(Dx)),Ex');
x = Q * xold;

%Pseudo-covariance 
pC = (x*transpose(x))/m;

% FIXED POINT ALGORITHM
 W = eye(n);
  Wold = zeros(n);
  k=0;
  
   while (norm(abs(Wold'*W)-eye(n),'fro')>(n*1e-5) && k < 15*n)
         k = k+1;
         Wold = W;
         
    for kk=1:n %Loop thru sources
        
        yy = W(:,kk)'*x;
       
        absy =abs(yy).^2; 
      %%Fixed point
      if type == 1 %%log
           g = 1./(a2 + absy);
           gp =  -1./(a2 + absy).^2;
      elseif type == 2  %Kurt
           g = absy;
           gp =  ones(size(absy));
      elseif type == 3  %sqrt
           g = 1./(2*sqrt(a2 + absy));
           gp =  -1./(4*(a2 + absy).^(3/2));
      end
       
       gRad =  mean(ones(n,1)*(g.*conj(yy)).*x,2);
       ggg = mean(gp.*absy + g);
       B = mean(gp .* conj(yy).^2)*pC;
       W(:,kk) =  Wold(:,kk)*ggg -(gRad) + (B*conj(Wold(:,kk)));
    
    end
    %Orthonormalization
    [E,D] = eig(W'*W);
    W = W * E * inv(sqrt(D)) * E';
   
  end; %Loop thru sources


shat = W'*x; %Estimated sources
Ahat = inv(Q)*W;
