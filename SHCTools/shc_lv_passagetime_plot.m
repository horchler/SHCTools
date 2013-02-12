function shc_lv_passagetime_plot(ti,varargin)
%SHC_LV_PASSAGETIME_PLOT  
%
%   SHC_LV_PASSAGETIME_PLOT(TI,TAU)
%   SHC_LV_PASSAGETIME_PLOT(TI,TP,TD)
%   SHC_LV_PASSAGETIME_PLOT(...,PLOTTYPE)
%   SHC_LV_PASSAGETIME_PLOT(...,C)

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-4-12
%   Revision: 1.0, 2-11-13


if nargin < 2
    error('SHCTools:shc_lv_passagetime_plot:TooFewInputs',...
          'Too few input arguments.');
elseif nargin == 2
    tau = varargin{1};
elseif nargin == 3
    tp = varargin{1};
    td = varargin{2};
else
    error('SHCTools:shc_lv_passagetime_plot:TooManyInputs',...
          'Too many input arguments.');
end

if iscell(ti)
    if ~all(cellfun(@(x)isvector(x) | isempty(x),ti))
        error('SHCTools:shc_lv_passagetime_plot:NonVectorCellTi',...
              '');
    end
    if ~all(cellfun(@isfloat,ti))
        error('SHCTools:shc_lv_passagetime_plot:NonFloatCellTi',...
              '');
    end
    
    if nargin == 2
        if ~all(cellfun(@(x)isvector(x) | isempty(x),tau))
            error('SHCTools:shc_lv_passagetime_plot:NonVectorCellTau',...
                  '');
        end
        if ~all(cellfun(@isfloat,tau))
            error('SHCTools:shc_lv_passagetime_plot:NonFloatCellTau',...
                  '');
        end
        if ~all(cellfun(@(x,y)isequal(numel(x),numel(y)),ti,tau))
            error('SHCTools:shc_lv_passagetime_plot:DimensionMismatchCellTau',...
                  '');
        end
        if all(cellfun(@(x)size(x,2) == 1 | isempty(x),tau))
            tau = vertcat(tau{:});
        else
            tau = horzcat(tau{:});
        end
        if ~isreal(tau) || ~all(isfinite(tau)) || ~all(tau > 0)
            error('SHCTools:shc_lv_passagetime_plot:InvalidCellTau',...
                  '');
        end
    else
        if ~all(cellfun(@(x)isvector(x) | isempty(x),tp))
            error('SHCTools:shc_lv_passagetime_plot:NonVectorCellTp',...
                  '');
        end
        if ~all(cellfun(@isfloat,tp))
            error('SHCTools:shc_lv_passagetime_plot:NonFloatCellTp',...
                  '');
        end
        if ~all(cellfun(@(x,y)isequal(numel(x),numel(y)),ti,tp))
            error('SHCTools:shc_lv_passagetime_plot:DimensionMismatchCellTp',...
                  '');
        end
        if all(cellfun(@(x)size(x,2) == 1 | isempty(x),tp))
            tp = vertcat(tp{:});
        else
            tp = horzcat(tp{:});
        end
        if ~isreal(tp) || ~all(isfinite(tp)) || ~all(tp > 0)
            error('SHCTools:shc_lv_passagetime_plot:InvalidCellTp',...
                  '');
        end
        
        if ~all(cellfun(@(x)isvector(x) | isempty(x),td))
            error('SHCTools:shc_lv_passagetime_plot:NonVectorCellTd',...
                  '');
        end
        if ~all(cellfun(@isfloat,td))
            error('SHCTools:shc_lv_passagetime_plot:NonFloatCellTd',...
                  '');
        end
        if ~all(cellfun(@(x,y)isequal(numel(x),numel(y)),ti,td))
            error('SHCTools:shc_lv_passagetime_plot:DimensionMismatchCellTd',...
                  '');
        end
        if all(cellfun(@(x)size(x,2) == 1 | isempty(x),td))
            td = vertcat(td{:});
        else
            td = horzcat(td{:});
        end
        if ~isreal(td) || ~all(isfinite(td)) || ~all(td > 0)
            error('SHCTools:shc_lv_passagetime_plot:InvalidCellTd',...
                  'td err');
        end
        tau = tp+td;
    end
    
    if length(ti) > 1
        nc = cumsum(cellfun(@numel,ti));
        multiColor = true;
    else
        multiColor = false;
    end
    if  all(cellfun(@(x)size(x,2) == 1 | isempty(x),ti))
        ti = vertcat(ti{:});
    else
        ti = horzcat(ti{:});
    end
    if ~isreal(ti) || ~all(isfinite(ti))
        error('SHCTools:shc_lv_passagetime_plot:InvalidCellTi',...
              '');
    end
    n = length(ti);
else
    if ~isvector(ti)
        error('SHCTools:shc_lv_passagetime_plot:NonVectorTi',...
              '');
    end
    if ~isreal(ti) || ~isfloat(ti) || ~all(isfinite(ti))
    	error('SHCTools:shc_lv_passagetime_plot:InvalidTi',...
              '');
    end
    n = length(ti);
    multiColor = false;
    
    if nargin == 2
        if ~isvector(tau)
            error('SHCTools:shc_lv_passagetime_plot:NonVectorTau',...
                  '');
        end
        if length(tau) ~= n
            error('SHCTools:shc_lv_passagetime_plot:DimensionMismatchTau',...
                  '');
        end
        if ~isreal(tau) || ~isfloat(tau) || ~all(isfinite(tau)) || ~all(tau > 0)
            error('SHCTools:shc_lv_passagetime_plot:InvalidTau',...
                  '');
        end
    else
        if ~isvector(tp)
            error('SHCTools:shc_lv_passagetime_plot:NonVectorTp',...
                  '');
        end
        if length(tp) ~= n
            error('SHCTools:shc_lv_passagetime_plot:DimensionMismatchTp',...
                  '');
        end
        if ~isreal(tp) || ~isfloat(tp) || ~all(isfinite(tp)) || ~all(tp > 0)
            error('SHCTools:shc_lv_passagetime_plot:InvalidTp',...
                  '');
        end
        
        if ~isvector(td)
            error('SHCTools:shc_lv_passagetime_plot:NonVectorTd',...
                  '');
        end
        if length(td) ~= n
            error('SHCTools:shc_lv_passagetime_plot:DimensionMismatchTd',...
                  '');
        end
        if ~isreal(td) || ~isfloat(td) || ~all(isfinite(td)) || ~all(td > 0)
            error('SHCTools:shc_lv_passagetime_plot:InvalidTd',...
                  '');
        end
        
        tau = tp+td;
    end
end

options = [];

% Open new figure if needed
cf = get(0,'CurrentFigure');
if isempty(cf) || strcmp(get(cf,'NextPlot'),'new')
    hf = figure('Renderer','zbuffer');
    ha = axes('Parent',hf);
else
    if isstruct(options) && isfield(options,'Parent')
        ha = options.Parent;
    elseif ishold
        ha = gca;
    else
        ha = newplot;
    end
    hf = ancestor(ha,'figure');
    set(hf,'Renderer','zbuffer');
end

% Size of figure axis in pixels
ax = get(ha,{'Children','Position','Units'});
p = hgconvertunits(hf,ax{2},ax{3},'pixels',hf);


% Half height of each frequency/period bar in pixels
h = 3;

% Adjust axes
if isempty(ax{1})
    mint = min(ti);
    [mt,mi] = max(ti);
    maxt = mt+tau(mi);
    maxfreq = 1./min(tau);
else
    xy = get(ha,{'XLim','YLim'});
    mint = min(min(ti),xy{1}(1));
    [mt,mi] = max(ti);
    maxt = max(mt+tau(mi),xy{1}(2));
    maxfreq = max(1./min(tau),xy{2}(2));
end

axis([mint maxt 0 maxfreq+h/p(4)])
tick = get(ha,{'XTick','YTick'});
dx = (tick{1}(end)-tick{1}(1))/(length(tick{1})-1);
dy = (tick{2}(end)-tick{2}(1))/(length(tick{2})-1);
mint(mint < tick{1}(1)) = tick{1}(1)-dx;
maxt(maxt > tick{1}(end)) = tick{1}(end)+dx;
maxfreq(maxfreq > tick{2}(end)) = tick{2}(end)+dy;
axis([mint maxt 0 maxfreq])

% Scaling to pixels to avoid round-off error in sizes and positions 
scx = p(3)/(maxt-mint);
scy = p(4)/maxfreq;

% Sort by time and save index
[tis,tii] = sort(ti);

% Colors
c = get(ha,'ColorOrder');
if multiColor
    for i = n:-1:1
        cp(1,i,:) = c(mod(find(tii(i) <= nc,1)-1,7)+1,:);
    end
else
    cp = c(1,:);
end

y = (1/scy)*bsxfun(@plus,round(scy./tau(tii(:,[1 1 1 1]))),[-h h h -h])';
rti = scx*tis;
if nargin == 3
    rtitp = rti+ceil(scx*tp(tii));
    rtitau = rtitp+ceil(scx*td(tii));
    
    x1 = (1/scx)*floor([[rti(1,[1 1]);rtitau(1:end-1,[1 1])] rtitp(:,[1 1])])';
    x2 = (1/scx)*floor([rtitp(:,[1 1]) rtitau(:,[1 1])])'; 
    patch(x1,y,cp,'EdgeColor','none')
    patch(x2,y,0.3*cp,'EdgeColor','none')
else
    rtitau = rti+ceil(scx*tau(tii));

    x1 = (1/scx)*floor([[rti(1,[1 1]);rtitau(1:end-1,[1 1])] rtitau(:,[1 1])])';
    patch(x1,y,cp,'EdgeColor','none')
end