function [U_multi]=multiinit_btd(d,d3,R,N,L,P)

NumRep   = 20; %Number of repetead initial analyses
options.MaxIter = 20; % Number of iterations in each initial fit
options.Display=1;
U=cpd_gevd(d3,R);
 for k=1:R
   tempimg=reshape(U{1,1}(:,k),60,[]);
   [U1,S,V] = svd(tempimg);
   U_gevd1{1,k}{1,1}=U1(:,1:L)*S(1:L,1:L);
   U_gevd1{1,k}{1,2}=V(:,1:L);
   U_gevd1{1,k}{1,4}=U{1,3}(:,k);
   U_gevd1{1,k}{1,3}=U{1,2}(:,k);
   U_gevd1{1,k}{1,5}=diag(ones(1,L(1)));
 end

[U,output] = btd_nls(d,U_gevd1,options); 
bestfit=output.fval(options.MaxIter+1);

for i=1:NumRep
    U1=lp_rnd(size(d),N,L,P);
    [U,output] = btd_nls(d,U1,options);
    if output.fval(options.MaxIter+1)<bestfit
        bestfit=output.fval(options.MaxIter+1);
        U_multi=U1;
    end
end
end