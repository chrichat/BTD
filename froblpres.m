function f = froblpres(T, U, varargin)
%FROBLPRES Frobenius norm of the residual for decomposition in lp terms
%   F = froblpres(T,U) and F = froblpres(T,U,N,L,P) compute the Frobenius 
%   norm of the residual tensor E which is given by lpgen(U)-T if U is 
%   given in the BTD format or lpgen(U,N,L,P)-T if U  is given in the CPD 
%   format.
%
%   froblpres(T,U,options), froblpres(T,U,N,L,P,options), froblpres(T,U,'key',
%   value) or froblpres(T,U,N,L,P,'key',value) can be used to set the following
%   options:
%
%   - ExpandLimit           Expands structured tensors to full tensors if the
%                           total number of entries prod(getsize(T)) is
%                           smaller than the value of this limit. Default:
%                           1e6.
%
%   All other options are passed on to lpres if needed.
%
%   See also lpres, froblvecres, lvecres, ll1res, frob, getstructure
    
% Author(s): Martijn Boussé      (Martijn.Bousse@esat.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/04/05   MB      Initial version

    if nargin >= 3 && ~ischar(varargin{1})
        N = varargin(1);
        varargin = varargin(2:end);
    else 
        N = {};
    end
    if nargin >= 4 && ~ischar(varargin{1})
        L = varargin(1);
        varargin = varargin(2:end);
    else 
        L = {};
    end 
    if nargin >= 5 && ~ischar(varargin{1})
        P = varargin(1);
        varargin = varargin(2:end);
    else 
        P = {};
    end    
    
    p = inputParser;
    p.addOptional('ExpandLimit', 1e6);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    passopt = p.Unmatched;
    
    type = getstructure(T);
    if any(strcmpi(type, {'full', 'incomplete', 'sparse'}))
        f = frob(lpres(T, U, N{:}, L{:}, P{:}, passopt));
    elseif prod(getsize(T)) <= options.ExpandLimit
        if ~isfield(passopt, 'Format'), passopt.Format = false; end
        f = frob(lpres(ful(T), U, N{:}, L{:}, P{:}, passopt));
    else 
        if all(cellfun(@isnumeric, U)) % CPD format
            if nargin < 4
                error('lpgen:parameters', ...
                      'N, L, and P should be given if U is in the CPD format.');
            end
            N = N{:}; L = L{:}; P = P{:};
            if any(cellfun('size',U,2) ~= [repmat(sum(L),1,N(1)), repmat(sum(P),1,N(2))])
                error('lp_nls:U', ...
                      ['cellfun(''size'',U,2) should be [sum(L), ... , sum(L),' ...
                       'sum(P), ..., sum(P).']);
            end 
            U = lpconvert(U,N,L,P);
        end
        f = frob(T,'squared') - 2*real(inprod(T,U)) + frob(U,'squared');
        f = sqrt(abs(f));
    end
end
