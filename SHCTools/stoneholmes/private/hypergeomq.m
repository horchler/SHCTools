function h=hypergeomq(n,d,z)
%HYPERGEOMQ  Fast numeric evaluation of the generalized hypergeometric function.
%   H = HYPERGEOMQ(N,D,Z) evaluates the generalized hypergeometric function for
%   the vector parameters N and D at the values in the array Z. N, D, and Z can
%   be any numeric or logical datatype and may be complex. The output, H, will
%   have the same dimensions as Z and will be double-precision unless Z is
%   single-precision, in which case, H will be as well.
%
%   HYPERGEOMQ uses the same low-level function as HYPERGEOM, but implements
%   several optimizations to achieve a performance boost of an order of
%   magnitude, or greater. Results from the two functions should agree to the
%   precision of the inputs. Additionaly, HYPERGEOMQ can avoid some errors due
%   to singularities that can occur when HYPERGEOM is evaluated with numeric Z.
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
%   See also: HYPERGEOM, SYM/HYPERGEOM.

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-22-12
%   Revision: 1.0, 8-11-12


% Check zeroes
if ~isvector(n) && ~isempty(n)
    error('SHCTools:hypergeomq:NonVectorN','The first input must be a vector.');
end
if ~isnumeric(n) && ~islogical(n)
    error('SHCTools:hypergeomq:NonNumericN',...
          'The first input argument must be a numeric or logical vector.');
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
if ~isnumeric(d) && ~islogical(d)
    error('SHCTools:hypergeomq:NonNumericD',...
          'The second input argument must be a numeric or logical vector.');
end
if ~all(isfinite(d))
    error('SHCTools:hypergeomq:NonFiniteD',...
          'The second input must be a vector of finite values.');
end

% Check Z values
if ~isnumeric(z) && ~islogical(z)
    error('SHCTools:hypergeomq:NonNumericD',...
          'The third input argument must be a numeric or logical array.');
end
if ~all(isfinite(z))
    error('SHCTools:hypergeomq:NonFiniteZ',...
          'The third input must be an array of finite values.');
end

if isa(z,'single')
    outType = 'single';
    scanStr = '%f32';
else
    outType = 'double';
    scanStr = '%f64';
end
if isempty(z)
    h = zeros(size(z),outType);
else
    % Convert zeroes to a string
    if isempty(n)
        N = '[]';
    else
        if isfloat(n)
            if isa(n,'single')
                printStr = '%.7g';
            else
                printStr = '%.15g';
            end
        else
            printStr = '%d';
        end
        if isreal(n)
            N = sprintf(printStr,n(1));
            if ~isscalar(n)
                N = ['[' N sprintf([',' printStr],n(2:end)) ']'];
            end
        else
            N = sprintf([printStr '+' printStr '*i'],real(n(1)),imag(n(1)));
            if ~isscalar(n)
                N = ['[' N sprintf([',' printStr '+' printStr '*i'],...
                    real(n(2:end)),imag(n(2:end))) ']'];
            end
        end
    end
    
    % Convert poles to a string
    if isempty(d)
        D = '[]';
    else
        if isfloat(d)
            if isa(d,'single')
                printStr = '%.7g';
            else
                printStr = '%.15g';
            end
        else
            printStr = '%d';
        end
        if isreal(d)
            D = sprintf(printStr,d(1));
            if ~isscalar(d)
                D = ['[' D sprintf([',' printStr],d(2:end)) ']'];
            end
        else
            D = sprintf([printStr '+' printStr '*i'],real(d(1)),imag(d(1)));
            if ~isscalar(d)
                D = ['[' D sprintf([',' printStr '+' printStr '*i'],...
                    real(d(2:end)),imag(d(2:end))) ']'];
            end
        end
    end
    
    % Set number of variable precision digits and format strings for writing z
    if isa(z,'double') || isa(z,'int64') || isa(z,'uint64')
        r = digits(17);
        printStr = '%.15f';
    elseif isa(z,'int32') || isa(z,'uint32')
        r = digits(10);
        printStr = '%.10f';
    else
        r = digits(9);
        printStr = '%.7f';
    end
    
    % Ensure the variable precision digits are reset to original value
    C = onCleanup(@()digits(r));
    
    % Try getting an analytic expression as a function of z
    mupadmexExists = (exist('mupadmex','file') == 3);
    
    if mupadmexExists
        % Low-level function with pre-formatted arguments, no validation
        sy = mupadmex('symobj::map','z','symobj::hypergeom',N,D);
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
        
        % Analytical expression can return errors, e.g. hypergeomq(1,0.5,-1)
        try
            h = eval(sc(2:end-1));
            return;
        catch	%#ok<CTCH>  
        end
    end
    
    % Convert Z values to a string
    if isreal(z)
        Z = sprintf(printStr,z(1));
        if ~isscalar(z)
            Z = ['[' Z sprintf([',' printStr],z(2:end)) ']'];
        end
    else
        Z = sprintf([printStr '+' printStr '*i'],real(z(1)),imag(z(1)));
        if ~isscalar(z)
            Z = ['[' Z sprintf([',' printStr '+' printStr '*i'],...
                real(z(2:end)),imag(z(2:end))) ']'];
        end
    end

    try
        if mupadmexExists
            % Low-level function with pre-formatted arguments and no validation
            hy = mupadmex('symobj::map',Z,'symobj::hypergeom',N,D);
            hc = mupadmex('symobj::char',charcmd(hy),0);
            if ~isempty(strncmp(hc,'"_symans_',9))
                hc = char(hy);
            end
        else
            % Pass string arguments and convert output to string
            hc = char(hypergeom(N,D,Z));
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
    
    % Replace special constants and format for conversion
    hc = strrep(hc,'RD_INF','Inf');
    hc = strrep(hc,'RD_NINF','-Inf');
    hc = strrep(hc,'RD_NAN','NaN');
    hc = strrep(hc,' ','');
    hc = strrep(hc,'*i','i');
    
    % Convert output to floating point - much faster than calling sym/double
    if isscalar(z)
        h = cast(str2double(hc(2:end-1)),outType);
    else
        % Trim "matrix([[ and ]])", format complex values as A + Bi
        hc = hc(11:end-4);
        if ~isempty(find(hc == 'i',1))
            hc = regexprep(hc,'([^,]*i)([^,]*)','$2+$1');
            hc = strrep(hc,'+-','-');
        end
        hc = textscan(hc,scanStr,'Delimiter',',');
        h = [hc{:}];

        % Reshape matrix and multi-dimensional array inputs
        if ~isvector(z)
            h = reshape(h,size(z));
        elseif ~isscalar(z) && size(z,1) == 1
            h = h.';
        end
    end
end