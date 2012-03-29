function adot=shc_lv_ode(t,a,varargin)	%#ok<INUSL>
%SHC_LV_ODE 
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 3-28-12
%   Revision: 1.0, 3-28-12


rho = varargin{1};
if nargin == 3
    mu = 0;
else
    mu = varargin{2};
end

a = min(max(a,0),1);
adot = a.*(rho(1:length(rho)+1:end)'-rho*a)+mu;