function T = lpblocks(N,L,P)
%LPCORE returns the core tensors of an lp decomposition in BTD format.

% Authors:  Martijn Boussé (Martijn.Bousse@esat.kuleuven.be)

    size_core = arrayfun(@(l,p) [l*ones(1,N(1)) p*ones(1,N(2))], L, P, 'UniformOutput', false);
    for r = 1:length(L)
        Tr.val = ones(1,L(r)*P(r));
        Tr.size = size_core{r};
        Tr.sub = [kron(ones(P(r),N(1)),(1:L(r)).'), kron((1:P(r)).',ones(L(r),N(2)))];
        Tr.sparse = 1;
        T{r} = ful(fmt(Tr));
    end
end

