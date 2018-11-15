function [hleg,hobj]=legendScatter(TextCell,SizeDim,Factor,LinSpec)
%% LEGENDSCATTER: workaround for R2016a and later for scatter plot legend.
%   This function aims to simplify the creation of a legend for scatter plots.
%   The function will plot a legend on current axis with the right
%   markersize proportion.
%   NB. From version 1.0.5 is added the compatibility with previous MATLAB
%       version.
%
%   Inputs:   - TextCell: a cell array containing the legend-entry text.
%               WARNING: THE FIRST CELL IS THE TITLE !!!
%             - SizeDim: the ''MarkerSize'' dimension for each entry.
%               WARNING: length(SizeDim)==length(TextCell)-1 !!!
%             - Factor: normaling factor for enlarge the marker size in
%                       a linear way. (MarkerSize*sqrt(Factor)
%             - LinSpec: specify the marker-string.
%
%   EXAMPLE: -----------------------------------------------------
%             Factor=1.5;
%             A=[1,1,1;
%                2,2,2;
%                3,3,3;
%                4,4,4;
%                5,5,5;
%                6,6,6;
%                7,7,7;
%                8,8,8;
%                9,9,9;
%                10,10,10];
% 
%             sh=scatter(A(:,1),A(:,2),(A(:,3).^2)*Factor,'ok');
%             %scattersize==MarkerSize*sqrt(Factor)pt)
%             leg=LEGENDSCATTER({'SIZE','1','2','3','4','5', ...
%                 '6','7','8','9','10'},[1:1:10],Factor,'ok');
%             leg.Location='NorthWest'; axis([0,11,0,11])
%            -----------------------------------------------------
%
%   VERSION:1.0.5
%   AUTHOR: Matteo Bagagli - ETH-Zurich // Oct. 2016
%   MAIL:   matteo.bagagli@erdw.ethz.ch

%% FAKE PLOTS
hold on
hplt=zeros(length(TextCell),1)';
for ii=1:length(TextCell)
    hplt(ii)=plot(NaN,NaN,LinSpec);
end
hold off

%% WORK
[hleg,hobj]=legend(hplt,TextCell,'Units','points');
drawnow % To update the figure and the legend
idx=(length(TextCell)+4):2:length(hobj); % Skip text and TITLE MARKER
v=version; 
if str2double(v(1)) >= 9                 % MATLAB v9.0 and higher
    % Title
    hobj(1).FontWeight='bold';
    ActPos=hobj(1).Position;
    hobj(1).Position=[0.5,ActPos(2),ActPos(3)];
    hobj(1).HorizontalAlignment='center';
    hobj(1).FontSize=11;
    hobj(length(TextCell)+2).Marker='none';
    % Body
    for ii=1:length(idx)
        hobj(idx(ii)).MarkerSize=SizeDim(ii)*sqrt(Factor);
    end
else                                     % MATLAB v8 and lower
    % Title
    set(hobj(1),'FontWeight','bold');
    ActPos=get(hobj(1),'Position');
    set(hobj(1),'Position',[0.5,ActPos(2),ActPos(3)]);
    set(hobj(1),'HorizontalAlignment','center');
    set(hobj(1),'FontSize',11);
    set(hobj(length(TextCell)+2),'Marker','none');
    % Body
    for ii=1:length(idx)
        set(hobj(idx(ii)),'MarkerSize',SizeDim(ii)*sqrt(Factor));
    end   
end % end switch version

end % EndMain
