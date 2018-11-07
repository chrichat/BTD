function T = lpgen(U, varargin)
%LPGEN Generate full tensor given as lp decomposition
%   T = lpgen(U) computes the tensor T given in the BTD format as the sum
%   of the R lp terms U{r} which are computed as tmprod(U{r}{N+1},
%   U{r}(1:N),1:N), where N is the order of the tensor T.
%
%   T = lpgen(U,N,L,P) computes the tensor T given in the CPD format as the
%   sum of the R lp terms by first converting to the BTD format.
%
%   T = lpgen(U,ind) for u in the BTD format and T = lpgen(U,N,L,P,ind) for
%   U in the CPD format compute only the elements of T corresponding to the
%   indices in ind. T has the same shape as the indices.
%
%   T = lpgen(U,i,j,k) for U in the BTD format and T = lpgen(U,L,i,j,k) for
%   U in the CPD format compute only the elements T(i,j,k) the tensor. The 
%   number of indices should match the order of the tensor given by 
%   length(U{1}) or length(U) for the BTD format or the CPD format, 
%   respectively. The colon operator can be used as a string.
%
%   See also cpdgen, ll1gen, lmlragen, btdgen, ttgen, lvecgen.

% Author(s): Martijn Boussé      (Martijn.Bousse@esat.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/09/03   MB      Initial version

if all(cellfun(@isnumeric, U))
    if nargin < 4
        error('lpgen:parameters', ...
              'N, L, and P should be given if U is in the CPD format.');
    end
    N = varargin{1}; 
    L = varargin{2}; 
    P = varargin{3};
    % ToDo: checks on N, L, and P.
    varargin = varargin(4:end);   
    T = btdgen(lpconvert(U,N,L,P),varargin{:});
else
    T = btdgen(U,varargin{:});
end


end

