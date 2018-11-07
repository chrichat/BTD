function E = lpres(T,U,varargin)
%LPRES Residual for a decomposition in lp terms
%   E = lpres(T,U) computes the residual tensor E as lpgen(U)-T, in which
%   U is given in the BTD format.
%
%   E = lpres(T,U,N,L,P) computes the residual tensor E as lpgen(U,N,L,P)-T, 
%   in which U is given in the CPD format.
%
%   lpres(T,U,options), lpres(T,U,N,L,P,options), lpres(T,U,'key',value) and
%   lpres(T,U,N,L,P,'key',value) can be used to set the following options:
%   
%     options.Format = [{true},false]   - If true, the tensor is formatted
%                                         using fmt(T) before computing the
%                                         residual.
%     options.Type                      - The type of the tensor. If not
%                                         given, the type is determined
%                                         automatically
%   See also lvecres, ll1res, cpdres, btdres, lmlrares
    
%   Authors: Martijn Boussé      (Martijn.Bousse@esat.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version History:
%   - 2016/04/04   MB      Initial version
  
    if all(cellfun(@isnumeric, U)) 
        if nargin < 5
            error('lpres:parameters', ...
                  'N, L, and P should be given if U is in the CPD format.');
        end
        N = varargin{1};
        L = varargin{2}; 
        P = varargin{3}; 
        % ToDo: checks on N, L, and P.
        varargin = varargin(4:end);
        
        if any(cellfun('size',U,2) ~= [repmat(sum(L),1,N(1)), repmat(sum(P),1,N(2))])
            error('lp_nls:U', ...
                  ['cellfun(''size'',U,2) should be [sum(L), ... , sum(L),' ...
                   'sum(P), ..., sum(P).']);
        end 
        
        E = btdres(T,lpconvert(U,N,L,P),varargin{:});
    else 
        E = btdres(T,U,varargin{:});
    end

end

