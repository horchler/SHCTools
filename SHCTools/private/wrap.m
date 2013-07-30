function p=wrap(q,cutoff,dim)
%WRAP  Wrap phase angles.
%   WRAP(Q) wraps radian phases Q by changing absolute jumps greater than or
%   equal to pi to their 2*pi complement. It wraps along the first non-singleton
%   dimension of Q. Q can be a scalar, vector, matrix, or N-D array.
%
%   WRAP(Q,TOL) uses a jump tolerance of TOL rather than the default TOL = pi.
%
%   WRAP(Q,[],DIM) wraps along dimension DIM using the default tolerance.
%   WRAP(Q,TOL,DIM) uses a jump tolerance of TOL.
%
%   Class support for input Q:
%       float: double, single
%
%   See also: UNWRAP, ANGLE, ABS, MOD

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-5-13
%   Revision: 1.0, 3-7-13


% Check Cutoff or set default
if nargin == 1 || isnumeric(cutoff) && isempty(cutoff)
    cutoff = pi;
else
    if ~isscalar(cutoff) || ~isnumeric(cutoff)
        error('wrap:InvalidCutoff','');
    end
end

if nargin < 3
    % Shift dimensions to work along first non-singleton dimension
    [q,nshifts] = shiftdim(q);
    
    % Wrap phase angles
    p = mod(q+cutoff,2*cutoff)-cutoff;
    
    % Shift back to original dimensions
    p = shiftdim(p,-nshifts);
else
    % Check Dimension
    if ~isscalar(dim) || ~isnumeric(dim)
        error('wrap:InvalidDimension','');
    end
    ndq = ndims(q);
    if dim < 1 || dim > ndq
        error('wrap:DimensionMismatch','');
    end
    
    % Shift dimensions to work along specified dimension
    perm = [dim:ndq 1:dim-1];
    q = permute(q,perm);
    
    % Wrap phase angles
    p = mod(q+cutoff,2*cutoff)-cutoff;
    
    % Shift back to original dimensions
    p = ipermute(p,perm);
end