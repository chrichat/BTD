function U = ll11_gevd_rank1(T,L,varargin)
%Ll11_GEVD  Initialization method for LL11 entirely based on LL1_GEVD

size_tens = size(T);

T_reshaped = reshape(T,[size_tens(1:2) prod(size_tens(3:4))]);
U_gevd_LL1 = ll1_gevd(T_reshaped,L,varargin{:});

U = cell(1,numel(U_gevd_LL1));
for r = 1:numel(U_gevd_LL1);
    
    % Rank 1 approximation
    CD = reshape(U_gevd_LL1{r}{3},size_tens(3:4));
    [CD_U,CD_S,CD_V] = svd(CD);
    
    a = CD_U(:,1)*CD_S(1,1);
    b = conj(CD_V(:,1));
    
    % Constructing new initialization
    U{r} = {U_gevd_LL1{r}{1},...
        U_gevd_LL1{r}{2},...
        a,...
        b,...
        U_gevd_LL1{r}{4}};
end

end