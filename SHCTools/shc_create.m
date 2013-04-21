function net=shc_create(nettype,params,varargin)
%SHC_CREATE  Create network and build RHO matrix from parameters.
%
%   NET = SHC_CREATE(NETTYPE,PARAMS)
%   NET = SHC_CREATE(NETTYPE,PARAMS,N)
%   NET = SHC_CREATE('custom',PARAMS,T)
%
%   See also:
%       SHC_CREATENETWORK, SHC_INITIALIZE, BUILDRHO, SHC_PARAMS, SHC_SYMPARAMS

%   Andrew D. Horchler, adh9@case.edu, Created 5-29-12
%   Revision: 1.0, 4-21-13


% Pass arguments directly to shc_createnetwork()
net = shc_createnetwork(nettype,params,varargin{:});
net = shc_initialize(net);
net = buildrho(net);