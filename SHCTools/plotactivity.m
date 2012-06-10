function plotactivity(a,t,net,opts)
%PLOTACTIVITY  Image plot of simulated SHC network activity.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 4-8-10
%   Revision: 1.0, 6-10-12


% Check activity data and time
if ~shc_ismatrix(a) || isempty(a) || ~isfloat(a)
    error('SHCTools:plotactivity:InvalidA',...
          'The activity data, A, must be a non-empty floating-point matrix.');
end
if ~isreal(a) || ~all(isfinite(a(:)))
    error('SHCTools:plotactivity:NonFiniteRealA',...
          'The activity data, A, must be a finite real floating-point matrix.');
end

if ~isvector(t) || isempty(t) || ~isfloat(t)
    error('SHCTools:plotactivity:InvalidT',...
          'The time, T, must be a non-empty floating-point vector.');
end
if ~isreal(t) || ~all(isfinite(t))
    error('SHCTools:plotactivity:NonFiniteRealT',...
          'The time, T, must be a finite real floating-point vector.');
end

if length(t) ~= size(a,1)
    error('SHCTools:plotactivity:DimensionMismatchT',...
         ['The length of the time vector, T, must be equal to the number of '...
          'rows of the activity matrix, A.']);
end

% Validate network structure
shc_validatenetwork(net);

n = size(a,2);
if n ~= net.size
    error('SHCTools:plotactivity:DimensionMismatchNetwork',...
         ['The number of columns of the the activity matrix, A, must be '...
          'equal to the ''size'' field of the SHC network.']);
end

y = 1:n;
yt = num2cell(sprintf('%d',y));
xtick = false;
cmap = jet(256);

% Handle options
if nargin > 3
    if ~iscell(opts) || ~(isempty(opts) && isnumeric(opts))
        if ischar(opts)
            opts = {opts};
        else
            error('SHCTools:plotactivity:NonCellOrStringOptions',...
                 ['Options argument must be a string or a cell arary of '...
                  'strings.']);
        end
    end
    for i = 1:length(opts)
        switch lower(opts{i})
            case 'nodenumbers'
                yt = 1:n;
            case 'timescale'
                xtick = true;
            case 'grayscale'
                cmap = 1-gray(256);
            otherwise
                error('SHCTools:plotactivity:UnknownOption',...
                      'Unknown option.');
        end
    end
end

yht = 0.5*y(end)/(n+1);
ymin = 0.5+yht/n;
ymax = y(end)+1-ymin;

% Set units and rescale time if needed
if t(end) <= 0.005
    tunits = '\mus';
    t = 1e6*t;
elseif t(end) <= 5
    tunits = 'ms';
    t = 1e3*t;
else
    tunits = 'sec.';
end

% Plot each state as an image
figure
hold on
for i = 1:n
    imagesc(t,y(i)+[1 -1]*yht,a(:,i)',[0 max(net.beta)])
end

% Set X-axis scaling and labeling
if xtick
    axis off
    set(gca,'XTickLabel',[])
    tickscale = 1.025;
    ymax = tickscale*ymax;
    
    % Draw time-scale bar
    taxis = 0.1*t(end);
    h = line([0 taxis],[0 0]+0.999*ymax);
    set(h,'Color',[0 0 0],'LineWidth',1.5)
    
    % Regex to remove insignificant trailing zeros from decimal values
    ts = regexprep(num2str(taxis,4),'(\..*[^0])0+|\.0+$()','$1');
    text(0,1.005*ymax,[ts ' ' tunits],'VerticalAlignment','top')
else
    set(gca,'Ycolor',[1 1 1],'YTick',y,'YTickLabel',yt)
    set(gca,'TickDir','out','TickLength',[0.05/(n+1) 0])
    xlabel(['Time (' tunits ')'])
end

% Y-label for each state ID
for i = 1:n
    text(-0.035*t(end),y(i),yt(i))
end
hy = ylabel('State ID');
set(hy,'Color',[0 0 0])

axis ij
axis([t(1) t(end) ymin ymax])
set(gcf,'Color',[1 1 1])
box off
colormap(cmap)

% Position colorbar
hc = colorbar;
if xtick
    p = get(gca,'Position');
    hcp = get(hc,'Position');
    set(hc,'Position',[hcp(1) p(2)+p(4)*(1-1/tickscale) hcp(3) p(4)/tickscale]);
    set(gca,'Position',p)
end
set(hc,'Box','off','Xcolor',[1 1 1],'TickDir','out')