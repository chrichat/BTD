   
size_tens = [10 10 10 10 10];
N = [3 2];
L = [2 3];
P = [1 2]; 

Ubtd = lp_rnd(size_tens,N,L,P);
Ucpd = lpconvert(Ubtd,N);

T = lpgen(Ubtd);
lpres(T,Ucpd,N,L,P);
froblpres(T,Ubtd,[1 1])

break

T(1,1,1) = NaN;
T = fmt(T);

options.MaxIter = 250;
options.Display = 1;
options.CGMaxIter = 50;
options.TolFun = eps^2;
options.TolX = eps;
% options.LargeScale = false;
[Uest,output] = lp_core(noisy(T,Inf),noisy(U,20),N,options);

% [Uest,output] = lp_core(T1,noisy(Ucpd,Inf),N,L,P);



