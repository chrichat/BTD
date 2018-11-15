function [U_multi]=multiinit(d,R)

NumRep   = 10; %Number of repetead initial analyses
options.MaxIter = 80; % Number of iterations in each initial fit


U_multi=cpd_gevd(d,R);
[U,output] = cpd_als(d,U_multi,options);
bestfit=output.fval(options.MaxIter+1);
for i=1:NumRep
    U1=cpd_rnd(size(d),R);
    [U,output] = cpd_als(d,U1,options);
    if output.fval(options.MaxIter+1)<bestfit
        bestfit=output.fval(options.MaxIter+1);
        U_multi=U1;
    end
end
end