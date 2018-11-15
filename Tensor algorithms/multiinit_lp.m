function [U_multi]=multiinit_lp(d,R,L,N)

NumRep   = 20; %Number of repetead initial analyses
options.MaxIter = 80; % Number of iterations in each initial fit
Lvar=(repmat(L,[1,R]));
P1=(repmat(1,[1,R]));
[ U_gevd1 ] = btd_gevd(d,R,L);
[U,output] = lp_nls(d,U_gevd1,N); 
bestfit=output.fval(options.MaxIter+1);
U_multi=U;
for i=1:NumRep
    U1=lp_rnd(size(d),N,Lvar,P1);
    [U,output] = lp_nls(d,U1,N,options);
    if output.fval(options.MaxIter+1)<bestfit
        bestfit=output.fval(options.MaxIter+1);
        U_multi=U1;
        sprintf('gevd is not the best')
    end
end
end