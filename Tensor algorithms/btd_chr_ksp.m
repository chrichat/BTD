function [ sol ] = btd_chr_ksp( data,L,R,varargin )
%   This function computed BTD as an expansion of multiple CPD with sdf
%   framework. L is the rank of BTD (different L for each mode has not yet
%   been implemented) and R the rank of the decomposition.
% When you give the optional parameter 'prior', it shall be followed by the knwon timecourse and a similarity 
% constant c 
% [ sol ] = btd_chr( data,L,R,'prior',di,c )
sizet=size(data);
options.Display = 0; % View convergence progress every 5 iterations.
options.TolFun=1e-9;
options.Tolx=eps;
nargin
if nargin<3
    error('The minimum number of inputs is 3 (data,L and R).')
end


for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i}),
      error (['Unknown type of optional parameter name (parameter' ...
	      ' names must be strings).']);
    end
    % change the value of parameter
    switch (varargin{i})
     case 'TolFun'
      options.TolFun = (varargin{i+1});
     case 'Display'
      options.Display = varargin{i+1};
     case 'Tolx'
      options.Tolx = varargin{i+1};
     case 'prior'
      prior=varargin{i+1};
      c=varargin{i+2};
      sizep=size(prior,2);
      %%%%DEFINE an anonymous function%%%
      fmri=@(z,task)struct_fmri(z,task,prior,c);
    end
  end
    for i=1:R
        model.variables.(sprintf('A%d',i)) =  randn(sizet(1),1);
        model.variables.(sprintf('B%d',i)) = randn(sizet(2),L);
        if exist('prior')
            if i>size(prior,2)
                model.variables.(sprintf('c%d',i))= randn(sizet(3),1); 
            else
                model.variables.(sprintf('prior%d',i))=prior(:,i);
            end
        else
            model.variables.(sprintf('c%d',i))= randn(sizet(3),1); 
        end
        if size(sizet,2)==4
            model.variables.(sprintf('d%d',i)) = randn(sizet(4),1); 
        end
    end

    model.factors.A1 = strsplit(sprintf('A%d ',ceil((1:L*R)/L)));
    model.factors.A1=model.factors.A1(1:end-1);
    model.factors.B1 = strsplit(sprintf('B%d ',1:R));
    model.factors.B1=model.factors.B1(1:end-1);
    if exist('prior')
       for i=1:L*sizep
           fprior{1,i}={sprintf('prior%d',ceil(i/L)),fmri};
       end       
       model.factors.C1=[fprior,strsplit(sprintf('c%d ',ceil(((L*sizep+1):L*R)/L)))];
       model.factors.C1=model.factors.C1(1:end-1);
    else
        model.factors.C1 = strsplit(sprintf('c%d ',ceil((1:L*R)/L)));
        model.factors.C1=model.factors.C1(1:end-1);
    end

    if size(sizet,2)==4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
         model.factors.D1 = strsplit(sprintf('d%d ',ceil((1:L*R)/L)));
         model.factors.D1=model.factors.D1(1:end-1);
    end
    
    model.factorizations.mybtd.data = data;
    if size(sizet,2)==4
        model.factorizations.mybtd.cpd = {'A1','B1','C1','D1'};
    else
        model.factorizations.mybtd.cpd = {'A1','B1','C1'};
    end
    % Solve the SDF model.
    
%     model.factorizations.myreg2.data = {zeros(sizet(2),L)};
%     model.factorizations.myreg2.regL1 = {'B1'};

    sol = sdf_nls(model,options);

end

