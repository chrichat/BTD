 function [U,output] = lp_rnd(size_tens,N,L,P,varargin)   
%LP_RND Pseudorandom initialization for lp decomposition.
%   U = lp_rnd(size_tens,N,L,P) generates R = length(L) = length(P)
%   pseudorandom terms U{r} to initalize algorithms that compute a lp
%   decomposition of an Nth-order tensor (N = N1+N2). Each term U{r} is a 
%   cell array of N1 factor matrices U{r}{n1} of sizes size_tens(n1) x L(r)
%   for n1 = 1, ..., N1, N2 factor matrices U{r}{n2} of sizes size_tens(n2)
%   x P(r) for n2 = 1, ..., N2, and an Nth-order core tensor U{r}{N+1} of
%   size L_{r} x ... x size L_{r} x size P_{r} x ... x size P_{r} with a
%   specific structure.
%
%   lp_rnd(T,N,L,P,) is shorthand for lp_rnd(getsize(T), N, L, P) if T is
%   real. If T is complex, then by default the terms will be generated
%   using pseudorandom complex numbers as well (cf. options).
%
%   lp_rnd(size_tens,N,L,P, options) and lp_rnd(T, L, P, N, options) may
%   be used to set the following options:
%
%      options.Real =         - The type of random number generator used to
%      [{@randn}|@rand|0]       generate the real part of the core tensor
%                               and factor matrices. If 0, there is no real
%                               part.
%      options.Imag =         - The type of random number generator used to
%      [@randn|@rand|0|...      generate the imaginary part of the core
%       {'auto'}]               tensor and factor matrices. If 0, there is
%                               no imaginary part. On 'auto', options.Imag
%                               is 0 unless the first argument is a complex
%                               tensor T, in which case it is equal to
%                               options.Real.
%      options.OutputFormat = - Format for the terms U: the BTD format
%                               returns R = length(L) = length(P)
%                               (rank-L(r) op rank-P(r)) terms with a
%                               specially structured core tensor, the CPD
%                               format returns N1 factor matrices with
%                               sum(L) columns and N2 factor matrices with
%                               sum(P) columns, respectively.
%
%   See also cpd_rnd, lmlra_rnd, btd_rnd, ll1_rnd, lvec_rnd
 
%   Authors: Martijn Boussé      (Martijn.Bousse@esat.kuleuven.be
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
 
    p = inputParser;
    p.addOptional('OutputFormat', 'btd');
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    passopt = p.Unmatched;

    if length(N) ~= 2
        error('lp_rnd:N','Length of N can only be 2.');
    end
    if sum(N) ~= length(size_tens)
        error('lp_rnd:N','Sum(N) should equal length(size_tens).');
    end
    if length(L) ~= length(P)
        error('lp_rnd:L','Length(L) should equal length(P).');
    end

    size_core = arrayfun(@(l,p) [l*ones(1,N(1)) p*ones(1,N(2))], L, P, 'UniformOutput', false);
    [U, output] = btd_rnd(size_tens, size_core, passopt);
    blocks = lpblocks(N,L,P);
    for r = 1:length(L), U{r}{end} = blocks{r}; end
    % Note: imposing specific structure here.

    if strcmpi(options.OutputFormat, 'cpd')
        U = lpconvert(U,N);
    elseif ~strcmpi(options.OutputFormat, 'btd')
        error('ll1_rnd:OutputFormat', 'Unknown output format %s', ...
              options.OutputFormat)
    end

end