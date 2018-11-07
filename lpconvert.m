function [Uout,varargout] = lpconvert(Uin,N,varargin)
%LPCONVERT convert decomposition in lp terms between CPD and BTD format.
%   UCP = lpconvert(UBTD,N) converts a decomposition in lp terms given in 
%   the BTD format to the CPD format. The BTD format UBTD is a cell with R 
%   = length(L) = length(P) terms, each having N1 factor matrices with L(r)
%   columns, N2 factor matrices with P(r) columns, and an Nth-order core
%   tensor of size L_{r} x ... x size L_{r} x size P_{r} x ... x size P_{r}
%   with a specific structure. We have that N = N1 + N2. The CPD format
%   UCPD is a cell with N factor matrices: N1 with sum(L) columns and N2
%   with sum(P) columns.
%
%   [UCPD,L,P] = lpconvert(UBTD,N) also returns the computed L and P.
%
%   UBTD = lpconvert(UCPD,N,L,P) converts a decomposition in lp terms
%   given in the CPD format to the BTD format.
%
%   lpconvert(U,N,'OutputFormat',format) or
%   lpconvert(U,N,L,P,'OutputFormat',format) can be used to explicitly set
%   the output format to 'btd' or 'cpd'.

%   Authors: Martijn Boussé (Martijn.Bousse@esat.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)


% ToDo: Derive N from the tensor instead of passing it as an argument?

    if nargin >= 3 && isnumeric(varargin{1})
        if ~isempty(varargin{1}), L = varargin{1}(:).'; end
        varargin = varargin(2:end);
    end
    if nargin >= 4 && isnumeric(varargin{1})
        if ~isempty(varargin{1}), P = varargin{1}(:).'; end
        varargin = varargin(2:end);
    end

    p = inputParser;
    p.addOptional('OutputFormat', 'auto');
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    Uin = Uin(:).';
       
    if all(cellfun(@isnumeric,Uin)) % CPD format
        if nargin < 4
            error('lpconvert:parameters', ...
                  'N, L, and P should be given if Ucpd is in the CPD format.');
        end
        if any(~cellfun(@ismatrix,Uin))
            error('lpconvert:Ucpd', ...
                  'Ucpd{n} should be matrices/vectors for all n.');
        end
        if any(cellfun('size',Uin,2) ~= [repmat(sum(L),1,N(1)), repmat(sum(P),1,N(2))])
            error('lpconvert:U', ...
                  ['cellfun(''size'',U,2) should be [sum(L), ... , sum(L),' ...
                   'sum(P), ..., sum(P).']);
        end
        
        inputFormat = 'cpd';
        if strcmpi(options.OutputFormat,'auto')
            options.OutputFormat = 'btd';
        end
    elseif all(cellfun(@iscell,Uin)) % BTD format
        Uin = cellfun(@(u) u(:).', Uin, 'UniformOutput', false);
        
        if length(N) ~= 2
            error('lp_rnd:N','Length of N can only be 2.');
        end
        size_tens = cellfun('size', Uin{1}(1:end-1), 1);
        if sum(N) ~= length(size_tens)
            error('lp_rnd:N','Sum(N) should equal the order of the tensor.');
        end

        Lest = cellfun(@(u) size(u{1},2), Uin);
        Pest = cellfun(@(u) size(u{N(1)+1},2), Uin);

        if exist('L','var')
            if length(Uin) ~= length(L)
                error('lpconvert:L', ['The detected number of terms length(Ubtd) ' ...
                                      'should match R = length(L).']);  
            end
            if any(L ~= Lest)
                error('lpconvert:L', ['The detected L and the given L are ' ...
                                      'not equal. (L is optional in the BTD format']);
            end
        else
            L = Lest;
            varargout{1} = L;
        end

        if exist('P','var')
            if length(Uin) ~= length(P)
                error('lpconvert:P', ['The detected number of terms length(Ubtd) ' ...
                                      'should match R = length(P).']);  
            end
            if any(P ~= Pest)
                error('lpconvert:P', ['The detected P and the given P are ' ...
                                      'not equal. (P is optional in the BTD format']);
            end
        else
            P = Pest;
            varargout{2} = P;
        end

        if any(cellfun(@(u) any(~cellfun(@ismatrix, u(1:end-1))), Uin))
            error('lpconvert:Ubtd', ['All entries Ubtd{r}{n}, n = 1, ..., N ' ...
                                       'should be matrices or vectors.']);
        end 
              
        if any(cellfun(@(u) size(u{1}, 2), cellfun(@(v) v(1:N(1)),Uin,'UniformOutput',false)) ~= L)
            error('lpconvert:Ubtd', ['All terms Ubtd{r}{n} should have L(r) ' ...
                                     ' columns, for n = 1:N(1).']);
        end
        if any(cellfun(@(u) size(u{1}, 2), cellfun(@(v) v(N(1)+1:sum(N)),Uin,'UniformOutput',false)) ~= P)
            error('lpconvert:Ubtd', ['All terms Ubtd{r}{n} should have P(r) ' ...
                                     ' columns, for n = N(1)+1:sum(N).']);
        end
        
        if any(arrayfun(@(v) any(cellfun(@(u) size(u{end},v),Uin) ~= L),1:N(1)))
            error('lpconvert:Ubtd', ['All terms U{r}{end} should be L(r) x ... ' ...
                                     'L(r) x P(r) x ... P(r) tensors of order N.']);
        end
        if any(arrayfun(@(v) any(cellfun(@(u) size(u{end},v),Uin) ~= P),N(1)+1:sum(N)))
            error('lpconvert:Ubtd', ['All terms U{r}{end} should be L(r) x ... ' ...
                                     'L(r) x P(r) x ... P(r) tensors of order N.']);
        end      
        if any(cellfun(@(u) any(cellfun('size',u(1:end-1),1)~=size_tens), Uin))
            error('lpconvert:Ubtd', ['All terms Ubtd{r} should define tensors ' ...
                                     'of the same size = cellfun(''size'',' ...
                                     'Ubtd{1}(1:end-1),1).']);
        end
            
        inputFormat = 'btd';
        if strcmpi(options.OutputFormat, 'auto')
            options.OutputFormat = 'cpd';
        end     
    else
        error('lpconvert:U','Unknown input format.');
    end
    
    if strcmpi(options.OutputFormat, 'btd')
        if strcmpi(inputFormat, 'cpd')
            % Do actual conversion
            U = Uin;
            Uin = [cellfun(@(u) mat2cell(u, size(u,1), L), U(1:N(1)),'UniformOutput',false) ...
                    cellfun(@(u) mat2cell(u, size(u,1), P), U(N(1)+1:sum(N)),'UniformOutput',false)];
            Uout = cell(1,length(L));
            blocks = lpblocks(N,L,P);
            for r = 1:length(L)
                Uout{r} = cellfun(@(u) u{r}, Uin, 'UniformOutput',false);
                Uout{r}{end+1} = blocks{r};
            end
        else
            Uout = Uin;
        end
        
    elseif strcmpi(options.OutputFormat, 'cpd')
        if strcmpi(inputFormat, 'btd')
            % Do actual conversion
            Uout = cell(1,length(Uin{1})-1);
            for n = 1:length(Uin{1})-1
                Uout{n} = cellfun(@(u) u{n}, Uin, 'UniformOutput', false);
                Uout{n} = cat(2, Uout{n}{:});
            end
        else
            Uout = Uin;
        end
        
    else
        error('lpconvert:OutputFormat', 'Unknown output format %s.', ...
                options.OutputFormat);
    end
    
    
end

