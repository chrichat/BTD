function [U,output] = lp_core(T,U0,N,varargin)
%LP_CORE Computational routines for lp decomposition.
%   lp_core should not be called directly. Use lp_minf or lp_nls instead.
%
%   [U,output] = lp_core(T,U0,N) decomposes the Nth-order tensor T in lp
%   terms starting from the initial solution U0. The result U is a cell of
%   length
%
%   [U,output] = lp_core(T,U0,N,P,L) can be used to specify an initial
%   solution in the CPD format. The solution U is also in the CPD format.
%
%   lp_core(T,U0,N,options), lp_core(T,U0,N,L,P,options),
%   lp_core(T,U0,N,'key',value) and lp_core(T,U0,N,L,P,'key',value) can be
%   used to set the following options:
%   
%      options.OutputFormat     - Either 'btd' or 'cpd'. If not given, the
%                                 same format as the input is used. 
%      options.OptimizationType - Either 'nls' (default) or 'minf' to select
%                                 second-order or first-order algorithms.
%      options.Algorithm        - The desired optimization method, e.g.
%                                 @nls_gndl, @nls_lm, @nls_gncg for NLS type
%                                 algorithms and @minf_lbfgsdl, @minf_lbfgs,
%                                 @minf_ncg for minf type algorithms. 
%      options.<...>            - Parameters passed to the selected method,
%                                 e.g., options.TolFun, options.TolX. See
%                                 also help [options.Algorithm].
%
%   See also: lp_minf, lp_nls

%   Authors: Martijn Boussé (Martijn.Bousse@esat.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

    type = getstructure(T);
    unstructuredtypes = {'full', 'incomplete', 'sparse'};
    isstructured = ~any(cellfun(@(s) strcmpi(type, s), ...
                                unstructuredtypes));
    if ~isstructured, 
        T = fmt(T,true); 
        type = getstructure(T);
    end
    size_tens = getsize(T);
    sz = getsize(U0);
    size_tens = [size_tens ones(1, length(sz)-length(size_tens))];
    
    isincomplete = strcmpi(type,'incomplete');

    if ~isempty(varargin) && isnumeric(varargin{1})
        L = varargin{1};
        varargin = varargin(2:end);
    else 
        L = [];
    end
    if ~isempty(varargin) && isnumeric(varargin{1})
        P = varargin{1};
        varargin = varargin(2:end);
    else 
        P = [];
    end
    L = L(:).';
    P = P(:).';
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('OutputFormat', 'auto');
    p.addOptional('LargeScale', 'auto');
    p.addOptional('Algorithm', @nls_gndl); % or @nls_gncgs,@nls_lm
    p.addOptional('OptimizationType', 'nls'); % or @nls_gncgs,@nls_lm
    p.addOptional('CGMaxIter', 10);
    p.addOptional('Display', 0);
    p.addOptional('TolLargeScale', 0.02);
    p.addOptional('TolAbs', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    
    % Set absolute tolerance
    cache = struct;
    cache.T2 = frob(T,'squared');
    if any(strcmpi(p.UsingDefaults, 'TolAbs')) 
        if isstructured
            options.TolAbs = 0.5*1e-15*cache.T2;
        else 
            options.TolAbs = 0.5*options.TolAbs*cache.T2;
        end
    end
    
    % convert to internal format
    if all(cellfun(@isnumeric, U0)) 
        % CPD format
        U0 = U0(:).';
        if isempty(L) || isempty(P) 
            error('lp_nls:parameters', ...
                  'L, and P should be given if U0 is in the CPD format');
        end
        if any(~cellfun(@ismatrix,U0))
            error('lp_nls:Ucpd', ...
                  'Ucpd{n} should be matrices/vectors for all n.');
        end
        if any(cellfun('size',U0,2) ~= [repmat(sum(L),1,N(1)), repmat(sum(P),1,N(2))])
            error('lp_nls:U', ...
                  ['cellfun(''size'',U,2) should be [sum(L), ... , sum(L),' ...
                   'sum(P), ..., sum(P).']);
        end       
        inputformat = 'cpd';
        R = length(L);
    elseif all(cellfun(@iscell, U0)) && ... 
            all(cellfun(@(u) all(cellfun(@isnumeric,u)), U0)) 
        % BTD format: convert to cpd format
        U0 = U0(:).';
        U0 = cellfun(@(u) u(:).', U0, 'UniformOutput', false);
        R = length(U0);
        [U0,L,P] = lpconvert(U0,N);    
        inputformat = 'btd';   
    else 
        error('lvec_core:U0', 'Unknown initialization format');
    end

    switch options.OutputFormat,
      case 'auto'
        options.OutputFormat = inputformat;
      case {'cpd', 'btd'}
        % ok
      otherwise
        error('lp_core:output_format', ['The output format should be ''auto'', ' ...
                            '''cpd'' or ''btd''']);
    end 
    
    sLP = sum(L.*P);
    
    if ischar(options.LargeScale)
        options.LargeScale = sum(size_tens)*sum(L.*P) > 100;
    end
    
    % Select solver
    options.OptimizationType = lower(options.OptimizationType);
    nlsfun = cellfun(@func2str, {@nls_gndl,@nls_gncgs,@nls_lm}, 'UniformOutput', ...
                     false);
    minffun = cellfun(@func2str, {@minf_lbfgsdl,@minf_lbfgs,@minf_ncg}, ...
                      'UniformOutput', false);
    switch options.OptimizationType
      case 'nls'
        if ~isfield(options, 'Algorithm')
            options.Algorithm = @nls_gndl;
        elseif ismember(func2str(options.Algorithm), minffun);
            error('lp:incompatibleAlgorithm', ...
                  ['The %s method is incompatible with the nls optimization ' ...
                   'type.'], func2str(options.Algorithm));
        elseif ~ismember(func2str(options.Algorithm), nlsfun)
            warning('lvec:unknownAlgorithm', ['The %s is not known by BTD. Use ' ...
                                'at your own risk'], options.Algorithm);
        end
      case 'minf'
        if ~isfield(options, 'Algorithm')
            options.Algorithm = @minf_lbfgsdl;
        elseif ismember(func2str(options.Algorithm), nlsfun)
            error('lp:incompatibleAlgorithm', ...
                  ['The %s method is incompatible with the minf optimization ' ...
                   'type.'], func2str(options.Algorithm));
        elseif ~ismember(func2str(options.Algorithm), minffun)
            warning('lp:unknownAlgorithm', ['The %s is not known by BTD. Use ' ...
                                'at your own risk'], options.Algorithm);
        end
      otherwise
        error('lp:unknownOptimizationType', ...
              ['The optimization type %s is unknown. Only ''nls'' and ''minf'' ' ...
               'are supported now'], options.OptimizationType);
    end
    
    % Select optimization subroutines.
    usestate = false;
    if strcmpi(options.OptimizationType, 'nls')
        usestate = true;
        if options.LargeScale, dF.JHJx = @JHJx; else dF.JHJ = @JHJ; end
        dF.JHF = @grad;
        if ~isfield(options, 'M'), options.M = 'block-Jacobi'; end
        switch options.M
          case 'block-Jacobi', dF.M = @M_blockJacobi;
          otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
        end
    else 
        dF = @grad;
    end
    state(U0,true);
    
    [U,output] = options.Algorithm(@objfun,dF,U0(:).',options);
    output.Name = func2str(options.Algorithm);
    
    if strcmpi(options.OutputFormat, 'btd')
        U = lpconvert(U,N,L,P,'OutputFormat','btd');
    end
    
    % For structured types, check the relative error
    if output.info == 4 && isstructured
        warning('lp_core:accuracy', ...
                ['Maximal numerical accuracy for structured tensors reached. The ' ...
                 'result may be improved using ful(T) instead of T.']);
    end
    
    function state(z,firstrun)
        
        if nargin == 2 && firstrun
            % Store the fraction of known elements.
            if isincomplete
                cache.scale = length(T.val)./prod(T.size);
            end
            cache.offset = cumsum([1 cellfun(@numel,z)]);
            expvec1 = arrayfun(@(n) repmat(sum(L(1:n-1))+(1:L(n)),1,P(n)), 1:R, 'UniformOutput', false);
            cache.expvec1 = cat(2, expvec1{:});
            expvec2 = arrayfun(@(n) kron(sum(P(1:n-1))+(1:P(n)),ones(1,L(n))), 1:R, 'UniformOutput', false);
            cache.expvec2 = cat(2, expvec2{:}); 
            P1 = arrayfun(@(n) repmat(eye(L(n)),P(n),1), 1:R, 'UniformOutput',false);
            cache.P1 = blkdiag(P1{:});
            P2 = eye(sum(P));
            cache.P2 = P2(:,cache.expvec2).';     
            cache.offset = cellfun(@numel, U0);
            cache.offset = cumsum([1 cache.offset]);
            if ~usestate, return; end
        end
        
        z = expand(z);
        cache.UHU = cellfun(@(u) conj(u'*u), z, 'UniformOutput', false);
        
        % Optionally cache some results for the block-Jacobi preconditioner.
        if ischar(options.M) || isa(options.M,'function_handle')
            UHU = cache.UHU;
            invUHU = cell(1,sum(N));
            for n = 1:N(1)
                idxw = true(1,sum(N));
                idxw(n) = false;
                invUHU{n} = contract1(inv(prod(cat(3,UHU{idxw}),3)),false);
            end
            for n = N(1)+1:sum(N)
                idxw = true(1,sum(N));
                idxw(n) = false;
                invUHU{n} = contract2(inv(prod(cat(3,UHU{idxw}),3)),false);
            end          
            cache.invUHU = invUHU;
        end   
    end

    function fval = objfun(z)
        % lp objective function.
        z = expand(z);
        if isstructured
            fval = abs(0.5*cache.T2 - real(inprod(T, z)) + 0.5*frob(z,'squared'));
        else 
            cache.residual = cpdres(T,z);
            if isstruct(cache.residual), fval = cache.residual.val;
            else fval = cache.residual; end
            fval = 0.5*(fval(:)'*fval(:)); %%%%% +l1(U1) norm %%%%
        end
    end

    function grad = grad(z)
        if usestate, state(z); end
        if ~isstructured, E = cache.residual; end
        offset = cache.offset;
        grad = nan(offset(end)-1, 1);
        z = expand(z);
        for n = 1:N(1)
            if isstructured,
                tmp = mtkrprod(z,z,n) - mtkrprod(T,z,n);
            else 
                tmp = full(mtkrprod(E, z, n)); %%% +sign(U1)
            end
            tmp = tmp*cache.P1;
            grad(offset(n):offset(n+1)-1) = tmp(:); 
        end 
        for n = N(1)+1:sum(N)
            if isstructured,
                tmp = mtkrprod(z,z,n) - mtkrprod(T,z,n);
            else 
                tmp = full(mtkrprod(E, z, n)); %%% +sign(U1)
            end
            tmp = tmp*cache.P2;
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end         
    end

    function jhj = JHJ(z)
        UHU = cellfun(@conj, cache.UHU, 'UniformOutput', false);
        offset = cache.offset;       
        jhj = zeros(offset(end)-1);
        z = expand(z);
        P1 = cache.P1;
        P2 = cache.P2;
        
        % Off diagonal blocks
        for n = 1:sum(N)
            idxw = true(1,sum(N));
            idxw(n) = false;
            idxn = offset(n):offset(n+1)-1;
            Wn = prod(cat(3, UHU{idxw}),3);            
            if n <= N(1)
                Wn = contract1(Wn,true);
            else
                Wn = contract2(Wn,true);
            end        
            jhj(idxn,idxn) = kron(Wn, eye(size_tens(n)));            
            for m = n+1:sum(N)
                idxw = true(1,sum(N));
                idxw([n m]) = false;
                idxn = offset(n):offset(n+1)-1;
                idxm = offset(m):offset(m+1)-1;
                Wnm = reshape(prod(cat(3,UHU{idxw}),3), [1 sLP 1 sLP]);
                JHJnm = bsxfun(@times, reshape(z{n}, [size_tens(n) 1 1 sLP]), ...
                               reshape(z{m}', [1 sLP size_tens(m) 1]));
                JHJnm = bsxfun(@times, JHJnm, Wnm);
                if n <= N(1) && m <= N(1), 
                    JHJnm = tmprod(JHJnm,{P1.',P1.'},[2 4]);
                    JHJnm = reshape(JHJnm, [size_tens(n)*sum(L) size_tens(m)*sum(L)]);
                elseif n <= N(1) && m > N(1), 
                    JHJnm = tmprod(JHJnm,{P1.',P2.'},[2 4]);
                    JHJnm = reshape(JHJnm, [size_tens(n)*sum(L) size_tens(m)*sum(P)]);
                elseif n > N(1) && m <= N(1),
                    JHJnm = tmprod(JHJnm,{P2.',P1.'},[2 4]);
                    JHJnm = reshape(JHJnm, [size_tens(n)*sum(P) size_tens(m)*sum(L)]);
                else
                    JHJnm = tmprod(JHJnm,{P2.',P2.'},[2 4]);
                    JHJnm = reshape(JHJnm, [size_tens(n)*sum(P) size_tens(m)*sum(P)]);
                end
                jhj(idxn,idxm) = JHJnm;
                jhj(idxm,idxn) = JHJnm';
            end 
        end
        
        % If incomplete, approximate the effect of missing entries.
        if isincomplete
            jhj = jhj*cache.scale;
        end
    end
    
    function y = JHJx(z,x)
        
        offset = cache.offset;
        UHU = cache.UHU;
        XHU = cell(1,sum(N));
        y = nan(offset(end)-1,1);

        % Diagonal blocks
        for n = 1:N(1)          
            idxw = true(1,sum(N));
            idxw(n) = false;
            Wn = prod(cat(3,UHU{idxw}),3); 
            Wn = contract1(Wn,true);          
            tmp = reshape(x(offset(n):offset(n+1)-1), size_tens(n), []);
            XHU{n} = conj(tmp'*z{n});
            y(offset(n):offset(n+1)-1) = tmp*Wn;           
        end
        for n = N(1)+1:sum(N)
            idxw = true(1,sum(N));
            idxw(n) = false;
            Wn = prod(cat(3,UHU{idxw}),3); 
            Wn = contract2(Wn,true);          
            tmp = reshape(x(offset(n):offset(n+1)-1), size_tens(n), []);
            XHU{n} = conj(tmp'*z{n});
            y(offset(n):offset(n+1)-1) = tmp*Wn; 
        end 
        XHUexp = expandfull(XHU);
                      
       % Off diagonal blocks       
        for n = 1:N(1) % upper block
            idxn = offset(n):offset(n+1)-1;
            Wn = zeros(sum(L.*P));
            for m = n+1:sum(N) 
                idxm = offset(m):offset(m+1)-1;
                idxw = true(1,sum(N)); 
                idxw([n m]) = false;                          
                Wn = Wn + prod(cat(3,UHU{idxw}),3).*XHUexp{m};
                Wnm = prod(cat(3,UHU{idxw}),3).*XHUexp{n};  
                if m <= N(1)
                    Wnm = contract1(Wnm,true);              
                else
                    Wnm = contract2(Wnm,true);
                end
                tmp = z{m}*Wnm;
                y(idxm) = y(idxm) + tmp(:);
            end
            Wn = contract1(Wn,true);
            tmp = z{n}*Wn;
            y(idxn) = y(idxn) + tmp(:);
        end      
        for n = N(1)+1:sum(N) % lower block
            idxn = offset(n):offset(n+1)-1;
            Wn = zeros(sum(L.*P));
            for m = n+1:sum(N) 
                idxm = offset(m):offset(m+1)-1;
                idxw = true(1,sum(N));
                idxw([n m]) = false;
                Wn = Wn + prod(cat(3,UHU{idxw}),3).*XHUexp{m};
                Wnm = prod(cat(3,UHU{idxw}),3).*XHUexp{n};
                if m <= N(1)
                    Wnm = contract1(Wnm,true);
                else
                    Wnm = contract2(Wnm,true);
                end
                tmp = z{m}*Wnm;
                y(idxm) = y(idxm) + tmp(:);
            end
            Wn = contract2(Wn,true);
            tmp = z{n}*Wn;
            y(idxn) = y(idxn) + tmp(:);
        end   
          
        % If incomplete, approximate the effect of missing entries.
        if isincomplete
            y = y*cache.scale;
        end
    end
    
    function x = M_blockJacobi(~,b)
    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
        x = nan(size(b));
        offset = cache.offset;
        invUHU = cache.invUHU;
        for n = 1:N(1)
            idx = offset(n):offset(n+1)-1; 
            x(idx) = reshape(b(idx), size_tens(n), sum(L))*invUHU{n};
        end
        for n = N(1)+1:sum(N)
            idx = offset(n):offset(n+1)-1; 
            x(idx) = reshape(b(idx), size_tens(n), sum(P))*invUHU{n};
        end
        
        % If incomplete, approximate the effect of missing entries.
        if isincomplete
            x = x/cache.scale;
        end
    end

    function W = contract1(W, transposed)
        P1 = cache.P1;       
        if transposed, W = P1.'*W*P1;
        else W = P1.'*W.'*P1; end
    end

    function W = contract2(W, transposed)
        P2 = cache.P2;       
        if transposed, W = P2.'*W*P2;
        else W = P2.'*W.'*P2; end
    end

    function V = expand(V)
        for n = 1:N(1)
            V{n} = V{n}(:, cache.expvec1);
        end
        for n = N(1)+1:sum(N)
            V{n} = V{n}(:, cache.expvec2);
        end
    end
 
    function V = expandfull(V)
       expvec1 = cache.expvec1;
       expvec2 = cache.expvec2;
       for n = 1:N(1)
            V{n} = V{n}(expvec1,expvec1);
        end
        for n = N(1)+1:sum(N)
            V{n} = V{n}(expvec2,expvec2);
        end
    end

end

