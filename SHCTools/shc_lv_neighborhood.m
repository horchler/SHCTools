function delta=shc_lv_neighborhood(net,delta_hat,epsilon,N)
%SHC_LV_NEIGHBORHOODSIZE    
%
%   DELTA = SHC_LV_NEIGHBORHOOD(NET,DELTA_HAT,EPSILON,N)
%
%   See also:
%       SHC_LV_INVPASSAGETIME

%   Andrew D. Horchler, adh9 @ case . edu, Created 6-1-12
%   Revision: 1.0, 2-12-13


% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_neighborhoodsize:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end
if ~shc_lv_iscycle(net)
    error('SHCTools:shc_lv_neighborhoodsize:NotCycle',...
          'The network must be an SHC cycle.');
end

t0 = 0;
dt = 1e-3;
tf = 1e3;
tspan = t0:dt:tf;

bet = net.beta(1);
opts = sdeset('EventsFUN',@(t,y)events(t,y,bet-delta_hat));

a0 = shc_lv_ic(net,delta_hat(1),epsilon);

[~,~,TE] = shc_lv_integrate(tspan,a0,net,epsilon,0,opts);
TE1 = diff(TE);
TE1 = TE1(1:2:end);

lte = length(TE1);
if lte < 3
    error('SHCTools:shc_lv_neighborhoodsize:TooFewCycles','Too few cycles');
end
M = ceil(N/lte);
TEj = cell(M,1);
TEj{1} = TE1;

parfor j = 2:M
    [~,~,TE] = shc_lv_integrate(tspan,a0,net,ep,0,opts);
    TE1 = diff(TE);
    TEj{j} = TE1(1:2:end);
end

delta = epsilon./shc_lv_invpassagetime(net,1,mean(vertcat(TEj{:})));



function [value,isterminal,direction]=events(t,y,d)	%#ok<INUSL>
value = y-d;
isterminal = 0;
direction = 0;