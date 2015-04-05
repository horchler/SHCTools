function rho=shc_lv_createcycle(alpha,beta,nu,direction)
%SHC_LV_CREATECYCLE  Create connection matrix RHO for Lotka-Volterra SHC cycle.
%   RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU) returns an N-dimensional
%   non-negative connection matrix RHO for the growth rates ALPHA, state
%   magnitudes BETA, and saddle values NU. ALPHA, BETA, and NU, must be symbolic
%   or floating-point scalars or length N vectors. If all three parameters are
%   scalar, a 3-by-3 connection matrix will be created. To form an SHC cycle,
%   the elements of RHO are:
%   
%                  { ALPHA(i)/BETA(i),                   if i == j
%       RHO(i,j) = { (ALPHA(i)-ALPHA(j)/NU(j))/BETA(j),  if i == mod(j,N)+1
%                  { (ALPHA(i)+ALPHA(j))/BETA(j),        otherwise
%   
%   where i,j are in {1,..., N}, N >= 3.
%   
%   RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU,-1) is the same as above except the
%   resultant connection matrix RHO is transposed, reversing the direction of
%   the cycle. RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU,1) is equivalent to
%   RHO = SHC_LV_CREATECYCLE(ALPHA,BETA,NU).
%   
%   See also:
%       SHC_LV_PARAMS, SHC_LV_EIGS, SHC_LV_JACOBIAN, SHC_LV_PASSAGETIME,
%       SHC_LV_MINTRANSITIONTIME

%   For details of the methods used, see:
%   
%   Andrew D. Horchler, Kathryn A. Daltorio, Hillel J. Chiel, and Roger D.
%   Quinn, "Designing Responsive Pattern Generators: Stable Heteroclinic Channel
%   Cycles for Modeling and Control," Bioinspiration & Biomimetics, Vol. 10,
%   No. 2., 2015, pp. 1-16.

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-27-14
%   Revision: 1.1, 4-5-15


% Dimension of parameters
n = max([length(alpha) length(beta) length(nu)]);
if n == 1
    n = 3;
end

% Expand scalars
z = ones(n,1);
alpha = alpha(:).*z;
beta = beta(:).*z;
nu = nu(:).*z;

% c >= 1, Adjust size of non-dominant stable eigenvalues, c == 1 for all equal
c = 1;

% Calculate Rho based on direction
if nargin > 3 && direction == -1
    rho = c*(alpha./beta)*z.';
    rho([n+1:n+1:end n]) = -alpha./(beta.*nu);
    rho(1:n+1:end) = 0;
    rho = rho+(z./beta)*alpha.';
else
    rho = c*z*(alpha./beta).';
    rho([2:n+1:end n*(n-1)+1]) = -alpha./(beta.*nu);
    rho(1:n+1:end) = 0;
    rho = rho+alpha*(z./beta).';
end

% Make non-negative
if ~isa(rho,'sym')
    rho(:) = max(rho(:),0);
end