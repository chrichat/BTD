 function [U,output] = cpdocpd_nls(T,U0,options)
%CPDOCPD_NLS CPDOCPD by nonlinear least squares.
%   [U,output] = cpdocpd_nls(T,U0) computes the factor matrices U{n}
%   corresponding to a decomposition that writes an Nth-order tensor T as a
%   sum of outer products of two CPDs by minimizing 
%   0.5*frob(T-cpdgen(U))^2; here, the factor matrices are expanded
%   accordingly. The algorithm is initialized with the terms matrices 
%   U0{r}. The structure output returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   cpdocpd_nls(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =   - The desired optimization method.
%      [@nls_gncgs| ...
%       {@nls_gndl}|@nls_lm]
%      options.M =           - The preconditioner to use when
%      [{'block-Jacobi'}|...   options.LargeScale is true.
%       false]
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX and
%                              options.PlaneSearchOptions. See also help
%                              [options.Algorithm].
%

%   Authors: Martijn Boussé (Martijn.Bousse@esat.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.
%   [3] ... ?
    
% Format the tensor T.
    T = fmt(T,true);
    if isstruct(T), size_tens = T.size;
    else size_tens = [size(T) ones(1,length(U0)-ndims(T))]; end % 
   
    R = length(U0);
    Nm = length(U0{1}{1});
    Ns = length(U0{1}{2});
    N = Nm+Ns;
    L = cellfun(@(u) size(u{1}{1},2), U0);
    P = cellfun(@(u) size(u{2}{1},2), U0);
    U = U0;
    for r = 1:R
        for nm = 1:Nm
            tmp = cellfun(@(u) u{1}{nm}, U, 'UniformOutput', false);
            U0{nm} = cat(2,tmp{:});
        end
    end
    for r = 1:R
        for ns = 1:Ns
            tmp = cellfun(@(u) u{2}{ns}, U, 'UniformOutput', false);
            U0{Nm+ns} = cat(2,tmp{:});
        end
    end
    
    %%
    expvecM = arrayfun(@(n) repmat(sum(L(1:n-1))+(1:L(n)),1,P(n)), 1:R, 'UniformOutput', false);
    expvecM = cat(2, expvecM{:});
    PM = arrayfun(@(n) repmat(eye(L(n)),P(n),1), 1:R, 'UniformOutput',false);
    PM = blkdiag(PM{:});

    expvecS = arrayfun(@(n) kron(sum(P(1:n-1))+(1:P(n)),ones(1,L(n))), 1:R, 'UniformOutput', false);
    expvecS = cat(2, expvecS{:});
    PS = eye(sum(P));
    PS = PS(:,expvecS).';

    % Check the options structure.
    isfunc = @(f)isa(f,'function_handle');
    xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
    if nargin < 3, options = struct; end
    if ~isfield(options,'Algorithm')
        funcs = {@nls_gndl,@nls_gncgs,@nls_lm};
        options.Algorithm = funcs{find(cellfun(xsfunc,funcs),1)};
    end
    if ~isfield(options,'CGMaxIter'), options.CGMaxIter = 10; end
    if ~isfield(options,'Display'), options.Display = 0; end
    if ~isfield(options,'M'), options.M = 'block-Jacobi'; end
    if ~isfield(options,'TolLargeScale'), options.TolLargeScale = 0.02; end

    % Call the optimization method.
    cache.offset = cellfun(@numel, U0);
    cache.offset = cumsum([1 cache.offset]);
    cache.T2 = frob(T,'squared');
    dF.JHJx = @JHJx;
    dF.JHF = @grad;
    switch options.M
      case 'block-Jacobi', dF.M = @M_blockJacobi;
      otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
    end
    state(U0,true);

    [U,output] = options.Algorithm(@objfun,dF,U0(:).',options);
    output.Name = func2str(options.Algorithm);
%     output.expvecM = expvecM;
%     output.expvecS = expvecS;
%     output.PM = PM;
%     output.PS = PS;
%     U = expand(U); % Other format
       
    function state(z,firstrun)       
        if nargin == 2 && firstrun
            % Store the fraction of known elements.
            if isstruct(T) && T.incomplete
                cache.scale = length(T.val)./prod(T.size);
            end
        end

        z = expand(z);
        UHU = zeros(sum(L.*P),sum(L.*P),N);
        for n = 1:N, UHU(:,:,n) = z{n}'*z{n}; end
        cache.UHU = UHU;
              
        % Optionally cache some results for the block-Jacobi preconditioner.
        if ischar(options.M) || isa(options.M,'function_handle')
            invUHU = cell(1,N);
            for n = 1:Nm
                invUHU{n} = contractM(inv(prod(UHU(:,:,[1:n-1 n+1:N]),3)), false);
%                 cond(invUHU{n})
            end
            for n = Nm+1:Nm+Ns
                invUHU{n} = contractS(inv(prod(UHU(:,:,[1:n-1 n+1:N]),3)), false);
%                 cond(invUHU{n})
            end
            cache.invUHU = invUHU;
        end
    end

    function fval = objfun(z)
    % cpdocpd objective function.
        z = expand(z);
        fval = cpdres(T,z);
        cache.residual = fval;
        fval = 0.5*(fval(:)'*fval(:));
    end

    function grad = grad(z)
        state(z);
        E = cache.residual;
        offset = cache.offset;
        grad = nan(offset(end)-1, 1);  
        z = expand(z);
        for n = 1:Nm
            tmp = mtkrprod(E, z, n);
            tmp = tmp*PM;
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
        for n = Nm+1:Nm+Ns
            tmp = mtkrprod(E, z, n);
            tmp = tmp*PS;
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
    end

    function y = JHJx(z,x)
        
        offset = cache.offset;
        UHU = cache.UHU;
        XHU = cell(1,N);
        y = nan(offset(end)-1,1);

        %% Diagonal blocks
        for n = 1:Nm
            Wn = prod(UHU(:,:,[1:n-1 n+1:N]),3);
            Wn = contractM(Wn);
            tmp = reshape(x(offset(n):offset(n+1)-1), [], sum(L));
            XHU{n} = conj(tmp'*z{n});
            y(offset(n):offset(n+1)-1) = tmp*conj(Wn);
        end
        for n = Nm+1:Nm+Ns
            Wn = prod(UHU(:,:,[1:n-1 n+1:N]),3);
            Wn = contractS(Wn);
            tmp = reshape(x(offset(n):offset(n+1)-1), [], sum(P));
            XHU{n} = conj(tmp'*z{n});
            y(offset(n):offset(n+1)-1) = tmp*conj(Wn);
        end   
                
        %% Off diagonal blocks
        ind = 1:N;
        for n = 1:Nm % left upper corner
            term = zeros(sum(L.*P),sum(L.*P));
            for m = [1:n-1 n+1:Nm]
                Wnm = prod(UHU(:,:,ind(~ismember(ind,ind([n m])))),3);
                term = term + conj(Wnm).*XHU{m}(expvecM, expvecM);
            end
            tmp = z{n}*contractM(term,true);
            idx = offset(n):offset(n+1)-1;
            y(idx) = y(idx) + tmp(:);
        end
        for n = 1:Nm % right upper corner
            term = zeros(sum(L.*P),sum(L.*P));
            for m = Nm+1:Nm+Ns
                Wnm = prod(UHU(:,:,ind(~ismember(ind,ind([n m])))),3);
                term = term + conj(Wnm).*XHU{m}(expvecS, expvecS);
            end
            tmp = z{n}*contractM(term,true);
            idx = offset(n):offset(n+1)-1;
            y(idx) = y(idx) + tmp(:);
        end  
        for n = Nm+1:Nm+Ns % left lower corner
            term = zeros(sum(L.*P),sum(L.*P));
            for m = 1:Nm
                Wnm = prod(UHU(:,:,ind(~ismember(ind,ind([n m])))),3);
                term = term + conj(Wnm).*XHU{m}(expvecM, expvecM);
            end
            tmp = z{n}*contractS(term,true);
            idx = offset(n):offset(n+1)-1;
            y(idx) = y(idx) + tmp(:);
        end   
        for n = Nm+1:Nm+Ns % right lower corner
            term = zeros(sum(L.*P),sum(L.*P));
            for m = [Nm+1:n-1 n+1:Nm+Ns]
                Wnm = prod(UHU(:,:,ind(~ismember(ind,ind([n m])))),3);
                term = term + conj(Wnm).*XHU{m}(expvecS, expvecS);
            end
            tmp = z{n}*contractS(term,true);
            idx = offset(n):offset(n+1)-1;
            y(idx) = y(idx) + tmp(:);
        end              
    end
    
    function x = M_blockJacobi(~,b)
    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
        x = nan(size(b));
        offset = cache.offset;
        invUHU = cache.invUHU;
        for n = 1:Nm
            idx = offset(n):offset(n+1)-1;
            x(idx) = reshape(b(idx), size_tens(n), sum(L))*invUHU{n};
        end
        for n = Nm+1:Nm+Ns
            idx = offset(n):offset(n+1)-1;
            x(idx) = reshape(b(idx), size_tens(n), sum(P))*invUHU{n};
        end        
    end
    
    function W = contractM(W, transposed)
        if nargin < 2, transposed = true; end
        if transposed
            W = PM.'*W*PM;
        else 
            W = PM.'*W.'*PM;
        end
    end

%     function W = contractMS(W, transposed)
%         if nargin < 2, transposed = true; end
%         if transposed
%             W = PM.'*W*PS;
%         else 
%             W = PM.'*W.'*PS;
%         end
%     end
% 
%     function W = contractSM(W, transposed)
%         if nargin < 2, transposed = true; end
%         if transposed
%             W = PS.'*W*PM;
%         else 
%             W = PS.'*W.'*PM;
%         end
%     end

    function W = contractS(W, transposed)
        if nargin < 2, transposed = true; end
        if transposed
            W = PS.'*W*PS;
        else 
            W = PS.'*W.'*PS;
        end
    end

     function V = expand(V)
        for n = 1:Nm
            V{n} = V{n}(:, expvecM);
        end
        for n = Nm+1:Nm+Ns
            V{n} = V{n}(:, expvecS);
        end
     end
 
%      function V = compress(V)
%         for n = 1:Nm
%             V{n} = V{n}*PM;
%         end
%         for n = Nm+1:Nm+Ns
%             V{n} = V{n}*PS;
%         end
%      end

end
