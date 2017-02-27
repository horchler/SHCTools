function h=hypergeomq(n,d,z)
%HYPERGEOMQ  Fast numeric evaluation of the generalized hypergeometric function.
%   H = HYPERGEOMQ(N,D,Z) evaluates the generalized hypergeometric function for
%   the vector parameters N and D at the values in the array Z. N, D, and Z can
%   be any numeric or logical datatype and may be complex. The output, H, will
%   have the same dimensions as Z and will be double-precision unless Z is
%   single-precision, in which case, H will be as well.
%
%   HYPERGEOMQ uses the same low-level function as HYPERGEOM, but implements
%   several optimizations to achieve a performance boost of approximately an
%   order of magnitude. Results from the two functions should agree to the
%   precision of the inputs. Additionaly, HYPERGEOMQ can avoid some errors due
%   to singularities that can occur when HYPERGEOM is evaluated with numeric Z
%   values.
%
%   Examples:
%       % An identity for the hypergeomentric function 1F0(a;;z)
%       z = -3:0.1:3; hz = hypergeomq(1,[],z); lz = 1./(1-z);
%       figure; plot(z,hz,'b',z,lz,'r.'); grid on; axis([-3 3 -10 10]);
%       xlabel('z'); ylabel('1F0(a;;z) == 1/(1-z)');
%
%       % Relationship between Confluent and Gaussian Hypergeometric Functions
%       z = logspace(1,2,20); h1 = hypergeomq(1,3,z); a2 = [2e2 5e2 2e3];
%       for i = 1:3; h2(i,:) = hypergeomq([1 a2(i)],3,z/a2(i)); end;
%       figure; loglog(z,abs(h1),'k',z,abs(h2)); grid on;
%       xlabel('z'); ylabel('|1F1(a;b;z)|, |2F1(a,b;c;z)|');
%       legend('_1F_1(a_1;b;z)','_2F_1(a_1,200;b;z/200)',...
%           '_2F_1(a_1,500;b;z/500)','_2F_1(a_1,2000;b;z/2000)',2);
%
%   Class support for ouput H:
%       float: double, single
%
%   See also: HYPERGEOM, SYM/HYPERGEOM.

%   Andrew D. Horchler, horchler @ gmail . com, Created 7-22-12
%   Revision: 1.2, 6-20-14


% Check zeroes
if ~isvector(n) && ~isempty(n)
    error('SHCTools:hypergeomq:NonVectorN','The first input must be a vector.');
end
if ~isnumeric(n) && ~islogical(n) && ~isa(n,'sym')
    error('SHCTools:hypergeomq:NonNumericN',...
         ['The first input argument must be a numeric, logical, or symbolic'...
          'vector.']);
end
if ~all(isfinite(n))
    error('SHCTools:hypergeomq:NonFiniteN',...
          'The first input must be a vector of finite values.');
end

% Check poles
if ~isvector(d) && ~isempty(d)
    error('SHCTools:hypergeomq:NonVectorD',...
          'The second input must be a vector.');
end
if ~isnumeric(d) && ~islogical(d) && ~isa(d,'sym')
    error('SHCTools:hypergeomq:NonNumericD',...
         ['The second input argument must be a numeric, logical, or symbolic'...
          'vector.']);
end
if ~all(isfinite(d))
    error('SHCTools:hypergeomq:NonFiniteD',...
          'The second input must be a vector of finite values.');
end

% Check Z values
if ~isnumeric(z) && ~islogical(z) && ~isa(z,'sym')
    error('SHCTools:hypergeomq:NonNumericZ',...
         ['The first input argument must be a numeric, logical, or symbolic'...
          'array.']);
end
if ~all(isfinite(z))
    error('SHCTools:hypergeomq:NonFiniteZ',...
          'The third input must be an array of finite values.');
end

isSym = isa(n,'sym') || isa(d,'sym') || isa(z,'sym');
if isSym
    dtype = 'sym';
elseif isa(z,'single')
    dtype = 'single';
else
    dtype = 'double';
end
if isempty(z)
    h = zeros(size(z),dtype);
else
    % Convert zeroes to a string
    if isa(n,'sym')
        N = sym2str(n);
    else
        N = float2str(n);
    end
    
    % Convert poles to a string
    if isa(d,'sym')
        D = sym2str(d);
    else
        D = float2str(d);
    end
    
    % Try getting an analytic expression as a function of z
    mupadmexExists = (exist('mupadmex','file') == 3);
    
    if mupadmexExists
        % Low-level function with pre-formatted arguments, no validation
        sy = mupadmex('symobj::hypergeom','z',N,D);
        sc = mupadmex('symobj::char',charcmd(sy),0);
        if ~isempty(strncmp(sc,'"_symans_',9))
            sc = char(sy);
        end
    else
        % Pass string arguments and convert output to string
        sc = char(hypergeom(N,D,'z'));
    end
    
    % Evaluate vectorized analytic expression if not function of hypergeometrics
    if isempty(strfind(sc,'hypergeom('))
        if ~isscalar(z)
            % Vectorize analytic expression
            sc = strrep(sc,'*','.*');
            sc = strrep(sc,'/','./');
            sc = strrep(sc,'^','.^');
        end
        
        % Analytical expression can return errors, e.g., hypergeomq(1,0.5,-1)
        try
            h = eval(sc);
            return;
        catch	%#ok<CTCH>  
        end
    end
    
    % Convert Z values to a string
    if isa(z,'sym')
        Z = sym2str(z);
    else
        Z = float2str(z);
    end

    % Evaluate numerically
    try
        if mupadmexExists
            % Low-level function with pre-formatted arguments and no validation
            h = mupadmex('map',Z,'symobj::hypergeom',N,D);
        else
            h = hypergeom(N,D,Z);
        end
    catch ME
        if strcmp(ME.identifier,'symbolic:mupadmex:CommandError')...
                && (~isempty(strfind(ME.message,'Singularity'))...
                || ~isempty(strfind(ME.message,'Division by zero')))
            error('SHCTools:hypergeomq:mupadmexSingularityError',...
                  'One or more singularities exist at or near the Z values.');
        else
            rethrow(ME);
        end
    end
    
    % Convert output to floating point - much faster than calling sym/double
    if ~isSym
        h = sym2float(h,dtype);
    end
    
    % Reshape matrix and multi-dimensional array inputs
    if ~isvector(z)
        h = reshape(h,size(z));
    elseif size(z,1) > 1
        h = h.';
    end
end