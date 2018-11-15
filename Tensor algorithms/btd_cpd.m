function [ sol ] = btd_cpd( T,R,L )
%%BTD based on sdf of 4 d tensor with cpd
size_tens=size(T)
for i=1:R
    eval(['model.variables.A' num2str(i) '= randn(size_tens(1),L)'])
    eval(['model.variables.B' num2str(i) '= randn(size_tens(2),L)'])
    eval(['model.variables.c' num2str(i) '= randn(size_tens(3),1)'])
    eval(['model.variables.d' num2str(i) '= randn(size_tens(4),1)'])
    temp=sprintf('A%d',i);
    model.factors.A{i}=temp;
    temp=sprintf('B%d',i);
    model.factors.B{i} =temp;
    temp=sprintf('c%d',i);
    temp1=sprintf('d%d',i);
    for j=1:L
        model.factors.C{(i-1)*L+j}=temp;
        model.factors.D{(i-1)*L+j}=temp1;
    end
end
model.factorizations.mybtd.data=T;
model.factorizations.mybtd.cpd={'A','B','C','D'};
options.Display=100;
options.TolFun=1e-12;
sol=sdf_nls(model,options);
end

