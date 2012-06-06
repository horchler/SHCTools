function net=shc_create(nettype,params,varargin)
%SHC_CREATE  Create network and build RHO matrix from parameters.
%
%   NET = SHC_CREATE(NETTYPE,PARAMS)
%   NET = SHC_CREATE(NETTYPE,PARAMS,N)
%   NET = SHC_CREATE('custom',PARAMS,T)

%   Andrew D. Horchler, adh9@case.edu, Created 5-29-12
%   Revision: 1.0, 5-29-12


% Pass arguments directly to shc_createnetwork()
net = shc_createnetwork(nettype,params,varargin{:});
net = shc_initialize(net);
net = buildrho(net);