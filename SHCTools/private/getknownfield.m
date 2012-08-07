function v=getknownfield(s,f,d)
%GETKNOWNFIELD	Get field f from struct s, or else yield default d.
%
%   See also:
%       SHC_LV_INTEGRATE

%   Andrew D. Horchler, adh9 @ case . edu, 10-28-10
%   Revision: 1.0, 6-30-12

if isfield(s,f)	% s could be empty.
	v = s.(f);
    if isempty(v)
        v = d;
    end
else
	v = d;
end