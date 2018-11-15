function [P,U0] = par_als(X,init,varargin)

%   P = par_als(X,R) computes an estimate of the best rank-R
%   CP model of a tensor X using an alternating least-squares
%   algorithm.  
%
%   P = par_als(X,U0) computes an estimate of the best length(U0)
%   CP model of a tensor X using an alternating least-squares
%   algorithm. 
%
%   P = par_als(X,...,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'toler' - Tolerance on difference in fit {1.0e-6}
%      'maxiters' - Maximum number of iterations {500}
%      'order' - Select the order through which the algorithm will loop
%      (in case a mode nees to remain same as the inital one dimension can
%      be obmitted)
%
%   [P,U0] = CP_ALS(...) also returns the initial guess.


%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = normtens(X);

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParamValue('toler',1e-6,@isscalar);
params.addParamValue('maxiters',500,@(x) isscalar(x) & x > 0);
params.addParamValue('dimorder',1:N,@(x) isequal(sort(x),1:N));
params.parse(varargin{:});

fittoler = params.Results.toler;
maxiters = params.Results.maxiters;
dimorder = params.Results.dimorder;

%% Compute order and/or Initial values
if iscell(init) ;  % Check the initial factor matrices U0.
    Uinit = init(:).'; 
    R = size(Uinit{1},2);
    if any(cellfun('size',Uinit,2) ~= R)
        error('cpd_als:U0','size(U0{n},2) should be the same for all n.');
    end
elseif rem(init,1) == 0 ; % then its an integer
    R=init;
    Uinit = cell(N,1);
      for n = 1:N
         Uinit{n} = rand(size(X,n),R);
      end
else
    error('The argument init is neither an order of number one terms nor initial matrices'); 
end

U = Uinit;
fit = 0;
B = zeros(R,R,N);
  for n = 1:N
     if ~isempty(U{n})
         B(:,:,n) = U{n}'*U{n};
      end
  end
    
for iter = 1:maxiters
   fitold = fit;
   % Iterate over all N modes of the tensor
   for n = dimorder(1:end)
        %create matriced form of tensor
        Xn = permute(X,[n 1:n-1,n+1:N]); 
        Xn = reshape(Xn, size(X,n), []);
        % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
        Z = khatrirao(U{[1:n-1,n+1:N]},'r');
        Unew = Xn*Z;
        
        V = prod(B(:,:,[1:n-1 n+1:N]),3);
        Unew = Unew / conj(V); 
        
        % Normalize each vector 
        if iter == 1
            lambda = sqrt(sum(Unew.^2,1))'; %2-norm
        else
            lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
        end              
        Unew = bsxfun(@rdivide, Unew, lambda');
        U{n} = Unew;
        B(:,:,n) = U{n}'*U{n};
   end
        
        P = reshape(U{1}*khatrirao(U(end:-1:2)).',cellfun('size',U(:).',1));
        U0=Uinit;
        residual=P-X;
        if normX == 0
            fit = normtens(P)^2 - 2 * innerprodtens(X,P);
        else
            normresidual = sqrt( normX^2 + normtens(P)^2 - 2 * innerprodtens(X,P) );
            fit = 1 - (normresidual / normX); %fraction explained by model
        end
        fitchange = abs(fitold - fit);
        
        % Check for convergence
        if (iter > 1) && (fitchange < fittoler)
            fprintf('Step size tolerance reached.');
            break;
        end
end   
fprintf('Maximum number of iterations reached.');
end

function n=normtens(T)
v = reshape(T, numel(T), 1);
n = norm(v);
end

function n=innerprodtens(X,Y)
   x = reshape(X, 1, numel(X));
   y = reshape(Y, numel(Y), 1);
   n = x*y;
end
   
