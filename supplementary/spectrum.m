function spectrum(A,name,varargin)
if nargin<2, name = 'A'; end

figure('name',['Spectrum ' name]);
if iscell(A)
    for i = 1:numel(A)
        plotspectrum(A{i},varargin{:});
        if i==1, hold all; end
    end
else
    A = squeeze(A);
    if size(A,3)>1
        l = cell(size(A,3),1);
        for i = 1:size(A,3)
            plotspectrum(A(:,:,i),varargin{:});
            if i==1, hold all; end
            l{i} = ['Slice ' num2str(i)];
        end
        legend(l)
    else
        plotspectrum(A,varargin{:});
    end
end

end

function plotspectrum(B,varargin)
svdB = svd(B);
if size(B,1)<size(B,2)
    svdB = [svdB;eps*ones(size(B,2)-size(B,1),1)];
end

if any(cellfun(@(x) strcmp(x,'Marker'),varargin)) || ...
        any(cellfun(@(x) strcmp(x,'LineStyle'),varargin))
    plot(svdB,varargin{:});
else
    plot(svdB,'*-',varargin{:});
end
xlabel('i'); ylabel('$\sigma_i$')
end