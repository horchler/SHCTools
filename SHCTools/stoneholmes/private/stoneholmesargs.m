function [szy,expnd]=stoneholmesargs(str,varargin)
%STONEHOLMESARGS  Checks input arguments for Stone-Holmes distribution functions
%   [SZY,EXPANSION] = STONEHOLMESARGS(STR,SZX,SZDELTA,SZEPSILON,SZLAMBDA_U)
%   returns the combined size vector, SZY, of the input size vectors and a
%   Boolean vector indicating if any of the input arrays are column vectors that
%   will need to be expanded to match SZY. STR is a string specifying the name
%   of the calling function in order to tailor error messages and check specific
%   constraints. Valid strings for this input format are 'pdf', 'cdf', 'inv',
%   'rnd', and 'like'.
%
%   [...] = STONEHOLMESARGS(STR,SZX,SZTHETA,SZLAMBDA_U) handles the two
%   parameter case, but only for the input strings 'pdf', 'cdf', and 'inv' as
%   'rnd' and 'like' do permit the parameters to be specified as Theta and
%   Lambda_U.
%
%   [...] = STONEHOLMESARGS(STR,SZDELTA,SZEPSILON,SZLAMBDA_U) and
%   [...] = STONEHOLMESARGS(STR,SZTHETA,SZLAMBDA_U) are as above, but valid
%   strings, STR, for this input format are 'median' and 'mode'.
%
%   [...] = STONEHOLMESARGS(STR,SZDELTA,SZEPSILON,SZLAMBDA_U,SZLAMBDA_S) and
%   [...] = STONEHOLMESARGS(STR,SZTHETA,SZLAMBDA_U,SZLAMBDA_S) are as above, but
%   the valid string, STR, for this input format is 'passagetime'.
%   
%   See also:
%       STONEHOLMESPDF, STONEHOLMESCDF, STONEHOLMESINV, STONEHOLMESRND,
%       STONEHOLMESLIKE, STONEHOLMESMEDIAN, STONEHOLMESMODE,
%       STONEHOLMESPASSAGETIME

%   ISROW, ISCOLUMN, and ISMATRIX are not used to maintain compatibility with
%   versions prior to Matlab 7.11 (R2010b).

%   Andrew D. Horchler, adh9 @ case . edu, Created 3-12-12
%   Revision: 1.0, 6-21-12


if ~ischar(str)
    error('SHCTools:stoneholmesargs:StringMissing',...
         ['First argument must be a string. Valid strings are: ''pdf'', '...
          '''cdf'', ''inv'', ''rnd'', ''like'', ''median'', ''mode'', and '...
          '''passagetime''.']);
end

% Special cases
isLike=strcmpi(str,'like');
isMedianMode=any(strcmpi(str,{'mode','median'}));
isLikeRnd=any(strcmpi(str,{'like','rnd'}));
isPassageTime=strcmpi(str,'passagetime');
isX=~(isMedianMode || isPassageTime);

if isMedianMode && nargin < 3 || isLikeRnd && nargin < 5 || nargin < 4
    error('SHCTools:stoneholmesargs:TooFewInputs','Not enough input arguments.')
end
if isMedianMode && nargin > 4 || nargin > 5
    error('SHCTools:stoneholmesargs:TooManyInputs','Too many input arguments.')
end

% Check for two parameter case
isTheta=true;
if isMedianMode && nargin == 4 || nargin == 5
    isTheta=false;
end

% Create size matrix used for comparisons
lv=cellfun('length',varargin);
if any(lv ~= lv(1))
    len=length(lv);
    sz=ones(len,max(lv));
    for i=1:len
        vi=varargin{i};
        sz(i,1:length(vi))=vi;
    end
else
    sz=cat(1,varargin{:});
end

% Non-scalar inputs
ns=~all(sz == 1,2);

% Check that first dimension is consistent for column expansion case
if ns(1) && any(sz([false;ns(2:end)],1) ~= sz(ns(1),1))
    if ~(isLike || isMedianMode)
        if isTheta
            thetastr='Theta'; 
        else
            thetastr='Delta, Epsilon';
        end
    end
    if isX && ns(1) && ~isLike
        if nargin < 5
            msgstr='the data';
        else
            switch(lower(str))
                case {'pdf','cdf'}
                    msgstr='X';
                case 'inv'
                    msgstr='P';
                case 'rnd'
                    msgstr='the specified dimensions';
                otherwise
                    msgstr='the data';
            end
        end
        errmsg=['If more than one of ' thetastr ', Lambda_U, or ' msgstr ...
                ' are non-scalar, they must have the same number of rows.'];
    elseif isPassageTime
        errmsg=['If more than one of ' thetastr ', Lambda_U, or Lambda_S '...
                'are non-scalar, they must have the same number of rows.'];
    elseif isLike
        if isTheta
            errmsg=['Theta and Lambda_U must be scalars or 1-by-N-by-... '...
                    'arrays.'];
        else
            errmsg=['Delta, Epsilon, and Lambda_U must be scalars or '...
                    '1-by-N-by-... arrays.'];
        end
    else
        if isTheta
            errmsg=['If Theta and/or Lambda_U are non-scalar, they must '...
                    'have the same number of rows.'];
        else
            errmsg=['If more than one of Delta, Epsilon, or Lambda_U, are '...
                    'non-scalar, they must have the same number of rows.'];
        end
    end
    error('SHCTools:stoneholmesargs:ParameterFirstDimensionMismatch',errmsg);
end

% Column vector and non-scalar, noncolumn vector array inputs
isCol=(sz(:,1) > 1 & all(sz(:,2:end) == 1,2));
isArray=(ns & ~isCol);

% Elements of size matrix that correspond to arrays
sza=sz(isArray,:);

% Check that all array (non-scalar, non-column vector) dimensions are consistent
if size(sza,1) > 1 && any(any(bsxfun(@ne,sza(2:end,2:end),sza(1,2:end)),2))
    if ~(isLike || isMedianMode)
        if isTheta
            thetastr='Theta'; 
        else
            thetastr='Delta, Epsilon';
        end
    end
    if isX && isArray(1)
        if nargin < 5
            msgstr='the data are non-scalar arrays and not column vectors';
        elseif isLike
            msgstr='X are non-scalar arrays';
        else
            switch(lower(str))
                case {'pdf','cdf'}
                    msgstr='X are non-scalar arrays and not column vectors';
                case 'inv'
                    msgstr='P are non-scalar arrays and not column vectors';
                case 'rnd'
                    msgstr=['the specified dimensions are non-scalar arrays '...
                            'and not column vectors'];
                otherwise
                    msgstr=['the data are non-scalar arrays and not column '...
                            'vectors'];
            end
        end
        errmsg=['If more than one of ' thetastr ', Lambda_U, or ' msgstr...
                ', they must have the same dimensions.'];
    elseif isPassageTime
        errmsg=['If more than one of ' thetastr ', Lambda_U, or Lambda_S '...
                'are non-scalar arrays and not column vectors, they must '...
                'have the same dimensions.'];
    elseif isLike
        if isTheta
            errmsg=['If Theta and/or Lambda_U are non-scalar arrays, they '...
                    'must have the same dimensions.'];
        else
            errmsg=['If more than one of Delta, Epsilon, or Lambda_U are '...
                    'non-scalar arrays, they must have the same dimensions.'];
        end
    else
        if isTheta
            errmsg=['If Theta and/or Lambda_U are non-scalar arrays and not '...
                    'column vectors, they must have the same dimensions.'];
        else
            errmsg=['If more than one of Delta, Epsilon, or Lambda_U are '...
                    'non-scalar arrays and not column vectors, they must '...
                    'have the same dimensions.'];
        end
    end
    error('SHCTools:stoneholmesargs:ParameterDimensionMismatch',errmsg)
end

% Set output size
if any(isArray)
    szy=sza(1,:);
elseif any(sz(:,1) == 0)
    szy=[0 1];
else
    if isLike
        szy=sz(1,:);
    else
        szy=max(sz);
    end
end
if isLike && szy(1) ~= 0
    szy(1)=sz(1);
end

% Set column vector expansion Boolean array
if nargout == 2
    if isLike
        expnd=(isCol(1) && any(isArray(2:end)));
    else
        if isX
            expnd=isCol & (1-eye(length(sz)))*isArray > 0;
            if nargin == 4
                expnd=[expnd(1);false;expnd(2:end)];
            end
        else
            if isPassageTime && nargin == 4 || isMedianMode && nargin == 3
                expnd=[false;(isCol & (1-eye(length(sz)))*double(isArray) > 0)];
            else
                expnd=(isCol & (1-eye(length(sz)))*double(isArray) > 0);
            end
        end
    end
end