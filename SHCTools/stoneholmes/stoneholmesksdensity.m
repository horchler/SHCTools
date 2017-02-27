function varargout=stoneholmesksdensity(tau,varargin)
%STONEHOLMESKSDENSITY  Compute kernel density for Stone-Holmes distribution.
%   [Y, TAUI] = STONEHOLMESKSDENSITY(TAU) computes a probability density
%   estimate of the sample in the vector TAU. The density estimate is evaluated
%   at at 100 points covering the range of the data. Y is the vector of density
%   values and TAUI is the set of 100 points. The estimate is based on a normal
%   kernel function, using a window parameter (bandwidth) that is a function of
%   the number of points in TAU. The Stone-Holmes distribution has support
%   on (0, Inf] so all values of TAU must be greater than zero.
%
%   Y = STONEHOLMESKSDENSITY(TAU,TAUI) specifies the vector TAUI of values where
%   the density estimate is to be evaluated. NaN is returned for any elements
%   of TAUI <= 0, i.e., outside of support of Stone-Holmes distribution.  
%
%   [Y, TAUI, U] = STONEHOLMESKSDENSITY(TAU,TAUI) also returns the bandwidth of
%   the kernel smoothing window.
%
%   [...] = STONEHOLMESKSDENSITY(...,'cdf') estimates the cumulative probability
%   density (CDF) instead of the PDF.
%
%   [...] = STONEHOLMESKSDENSITY(...,'icdf') estimates the inverse cumulative
%   probability density instead of the PDF. For [...] = ksdensity(TAU,YI,'icdf')
%   the inverse CDF of the values in TAU is estimated and evaluated at the
%   probability values specified in YI.
%
%   Example:
%       % Compare density of random data to exact Stone-Holmes distribution PDF
%       delta=1; epsilon=0.01; lambda_u=0.5; lambda_s=1; n=1e3;
%       r=stoneholmesrnd(delta,epsilon,lambda_u,lambda_s,n,1);
%       x1=0:0.1:25; binx=0.5; edges=x1(1):binx:x1(end);
%       n=histc(r,edges-0.5*binx); p=n/(binx*sum(n));
%       y1=stoneholmespdf(x1,delta,epsilon,lambda_u,lambda_s);
%       [y2,x2]=stoneholmesksdensity(r,x1);
%       figure; bar(edges,p,1); hold on; plot(x1,y1,'g',x2,y2,'c');
%       shading flat; h=xlabel('$x$'); h(2)=ylabel('p($x$)');
%       h(3)=title('Kernel Density for Stone-Holmes distribution');
%       set(h,'Interpreter','latex'); axis([0 25 0 0.3]);
%       
%   See also:
%       KSDENSITY, STONEHOLMESKSTEST, STONEHOLMESCHI2GOF, STONEHOLMESFIT,
%       STONEHOLMESCDF, STONEHOLMESPDF, STONEHOLMESRND, STONEHOLMESLIKE,
%       STONEHOLMESPASSAGETIME

%   Andrew D. Horchler, horchler @ gmail . com, Created 6-5-13
%   Revision: 1.1, 6-5-13


if nargout > 4
    error('SHCTools:stoneholmesksdensity:TooManyOutputs',...
          'Too many output arguments.');
end

% Check Tau
if ~isvector(tau) || ~isfloat(tau) || any(tau <= 0)
    error('SHCTools:stoneholmesksdensity:InvalidTau',...
          'TAU must be a floating point vector of strictly positive values.');
end

% Check variable inputs and compute kernel density
if nargin > 3
    error('SHCTools:stoneholmesksdensity:TooManyInputs',...
          'Too many input arguments.');
elseif nargin == 3
    taui = varargin{1};
    if ~isvector(taui) || ~isfloat(taui) || isscalar(taui) && taui <= 0
        error('SHCTools:stoneholmesksdensity:InvalidPoints2',...
             ['TAUI must be a floating point vector. (If TAUI is a scalar, '...
              'then it must be greater than zero.']);
    end
    % Set TAUI elements outside support of Stone-Holmes distribution to NaN
    taui(taui <= 0) = NaN;
    
    fun = varargin{2};
    if ~ischar(fun) || ~any(strcmp(fun,{'cdf','icdf'}))
        error('SHCTools:stoneholmesksdensity:InvalidFunction2',...
             ['The optional function type must be ''pdf'' (default), '...
              '''cdf'', or ''icdf''.']);
    end
    
    % Compute kernel density
    [varargout{1:min(max(nargout,1),4)}] = ksdensity(tau,taui,'function',fun,...
                                                     'support','positive');
elseif nargin == 2
    if ischar(varargin{1})
        fun = varargin{1};
        if ~any(strcmp(fun,{'cdf','icdf','pdf'}))
            error('SHCTools:stoneholmesksdensity:InvalidFunction1',...
                 ['The optional function type must be ''pdf'' (default), '...
                  '''cdf'', or ''icdf''.']);
        end
        
        % Compute kernel density
        [varargout{1:min(max(nargout,1),4)}] = ksdensity(tau,'function',fun,...
                                                         'support','positive');
    else
        taui = varargin{1};
        if ~isvector(taui) || isscalar(taui) && taui <= 0
            error('SHCTools:stoneholmesksdensity:InvalidPoints1',...
                 ['TAUI must be a floating point vector. (If TAUI is a '...
                  'scalar, then it must be greater than zero.']);
        end
        % Set TAUI elements outside support of Stone-Holmes distribution to NaN
        taui(taui <= 0) = NaN;
        
        % Compute kernel density
        [varargout{1:min(max(nargout,1),4)}] = ksdensity(tau,taui,...
                                                         'support','positive');
    end
else
    % Compute kernel density
    [varargout{1:min(max(nargout,1),4)}] = ksdensity(tau,'support','positive');
end