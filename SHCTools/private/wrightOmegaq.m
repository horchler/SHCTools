function W=wrightOmegaq(Z)
%WRIGHTOMEGAQ  Wright omega function, a solution of the equation W+LOG(W) = Z.
%   W = WRIGHTOMEGAQ(Z) performs floating point evaluation of the Wright omega
%   function. Z is an array and may be complex. If Z is an array of symbolic
%   values, it is converted to double-precision for computation and then recast
%   as symbolic.
%
%   Example:
%       % Plot magnitude of the function surface over the complex plane
%     	v = -6*pi:3*pi/25:6*pi; [x,y] = meshgrid(v);
%       w = wrightOmegaq(x+y*1i);
%
%       figure('Renderer','zbuffer'); surf(v,v,abs(w)); colormap(hot(256));
%       shading flat; axis square; xlabel('X'); ylabel('Y');
%       zlabel('|\omega(X+Yi)|'); title('Wright \omega Function')
%
%   Note:
%       Due to numerical precision and the nature of this equation, the inverse
%       of the Wright omega funtion, may not always return a value close (in an
%       absolute sense) to Z, e.g., Z <= -714.84989998498+Y*1i, -pi < Y <= pi.
%
%   Class support for Z:
%       float: double, single
%       symbolic
%
%   See also: WRIGHTOMEGA, LAMBERTW

%   Based on:
%
%   [1] Piers W. Lawrence, Robert M. Corless, and David J. Jeffrey, "Algorithm
%   917: Complex Double-Precision Evaluation of the Wright omega Function," ACM
%   Transactions on Mathematical Software, Vol. 38, No. 3, Article 20, pp. 1-17,
%   Apr. 2012. http://dx.doi.org/10.1145/2168773.2168779
%
%	[2] Robert M. Corless and David J. Jeffrey, "The Wright omega Function," In:
%   Artificial Intelligence, Automated Reasoning, and Symbolic Computation,
%   Joint International Conferences, AISC 2002 and Calculemus 2002, Marseille,
%   France, July 2002, (Jacques Calmet, Belaid Benhamou, Olga Caprotti, Laurent
%   Henocque, and Volker Sorge, Eds.), Berlin: Springer-Verlag, pp. 76-89, 2002.
%   http://orcca.on.ca/TechReports/2000/TR-00-12.html
%
%	Numbers in parentheses below refer to equations in Lawrence, et al. 2012.

% 	The inverse Wright omega function is defined as (Corless & Jeffrey, 2002):
%       Z = W+LOG(W)-2*pi*1i,  	-Inf < REAL(W) < -1, IMAG(W) = 0
%       Z = -1+/-pi*1i,         W = -1
%       Z = W+LOG(W),        	otherwise

%   The solution of W+LOG(W) = Z is given by (Corless & Jeffrey, 2002):
%       W = WRIGHTOMEGA(Z),                      	Z ~= t+/-pi*1i, t <= -1
%       W = WRIGHTOMEGA(Z),WRIGHTOMEGA(Z-2*pi*1i),	Z = t+pi*1i, t <= -1
%       W = NaN,                                  	Z = t-pi*1i, t <= -1

%   WRIGHTOMEGAQ is up three to four orders of magnitude faster than WRIGHTOMEGA
%   for double-precision arrays. Additionally, it has much less numeric error,
%   properly evaluates values greater than 2^28, supports single-precision
%   evaluation, and handles NaN inputs.

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-12-12
%   Revision: 1.0, 3-12-12


% Convert symbolic input, converted back at end
isSym = isa(Z,'sym') || ischar(Z);
if isSym
    try
        Z = double(sym(Z));
    catch ME
        if strcmp(ME.identifier,'symbolic:sym:double:cantconvert')
            error('SHCTools:wrightOmegaq:InvalidSymbolicZ',...
                 ['Symbolic input Z must contain only numeric values and no '...
                  'expressions containing variables.'])
        else
            rethrow(ME);
        end
    end
elseif ~isfloat(Z)
    error('SHCTools:wrightOmegaq:InvalidZ',...
          'Input Z must be an array of floating point or symbolic values.');
end

% Support for single precision: single(pi) ~= pi
dataType = class(Z);
PI = cast(pi,dataType);
tol = eps(dataType);

X = real(Z(:));
Y = imag(Z(:));

if isempty(Z) || all(isnan(Z(:)))
    W = Z;
elseif isscalar(Z)
    % Special values
    if Z > 2^59
        W = Z;                  % W self-saturates: X > 2^59 (abs(Y) > 2^54 too)
    elseif Z == 0
        W(1) = 0.5671432904097838;                              % Omega constant
    elseif Z == 1
        W(1) = 1;
    elseif Z == 1+exp(1)
        W(1) = exp(1);
    elseif isreal(Z) || Y == 0
        if Z < log(eps(realmin(dataType)))-log(2)
            W(1) = 0;                                               % Z -> -Inf
        else
            % W(1) used in order retain datatype
            if Z <= -2
                % Region 3: series about -Inf
                x = exp(X);
                W(1) = x*(1-x*(1-x*(36-x*(64-125*x))/24));          	% (24)
                
                % Series is exact, X < -exp(2)
              	if X < -7.38905609893065
                    return;
                end
            elseif Z > pi+1
                % Region 7: log series about Z = Inf
                x = log(X);
                lzi = x/X;
                W(1) = X-x+lzi*(1+lzi*(x/2-1+lzi*((x/3-3/2)*x+1)));   	% (25)
            else
                % Region 4: series about Z = 1
                x = X-1;
                W(1) = 1+x*(1/2+x*(1/16-x*(1/192+x*(1/3072 ...
                       -(13/61440)*x))));                           	% (29)
            end
            
            % Residual
            r = X-(W+log(W));                                           % (14)
            
            if abs(r) > tol
                % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
                W = W*(1+(r/(1+W))*((1+W)*(1+W+(2/3)*r)...
                    -r/2)/((1+W)*(1+W+(2/3)*r)-r));                     % (15)

                % Test residual
                r = X-(W+log(W));                                       % (14)

                % Second iterative improvement via FSC method, if needed
                if abs(r) > tol
                    W = W*(1+(r/(1+W))*((1+W)*(1+W+(2/3)*r)...
                        -r/2)/((1+W)*(1+W+(2/3)*r)-r));             	% (15)
                end
            end
        end
    else
        % Special complex values
        if Z < log(eps(realmin(dataType)))-log(2) && Y > -PI && Y <= PI
            W(1) = 0;                                               % Z -> -Inf
        elseif Z == -1+PI*1i || Z == -1-PI*1i
            W(1) = -1;
        elseif Z == log(1/3)-1/3+PI*1i
            W(1) = -1/3;
        elseif Z == log(2)-2-PI*1i
            W(1) = -2;
        elseif Z == (1+PI/2)*1i
            W(1) = 1i;
        else
            xgtm2 = (X > -2);
            
            % W(1) used in order retain datatype
            if X <= -1 && abs(Y) == PI
                % Upper and lower lines of discontinuity
                if Y == PI
                    % Region 3: series about -Inf
                    x = exp(X);
                    W(1) = ((((-125*x-64)*x-36)*x/24-1)*x-1)*x;      	% (24)

                    % Series is exact, X < -exp(2)
                    if X < -7.38905609893065
                        return;
                    end
                elseif xgtm2
                    % Regions 1 and 2: near Z = -1+pi*1i and Z = -1-pi*1i
                    x = sign(Y)*sqrt(-2*(X+1));                   	% (20, 22)
                    W(1) = -1+x*(1-x*(1440-x*(120+x*(16+x)))/4320);	% (21, 23)
                else
                    % Region 6: negative log series about Z = Z+pi*1i
                    x = log(-X);
                    lzi = x/X;
                    W(1) = X-x+lzi*(1+lzi*(x/2-1+lzi*((x/3-3/2)*x+1)));	% (28)
                end
                
                % Regularization: adjust Z and flip sign of W
                Z = X;
                s = -1;
            else
                xlteq1 = (X <= 1);
                r12x = (xgtm2 && xlteq1);
                
                if ~xgtm2 && Y > -PI && Y < PI
                    % Region 3: series about -Inf, within lines of discontinuity
                    x = exp(Z);
                    W(1) = x*(1-x*(1-x*(36-x*(64-125*x))/24));        	% (24)
                    
                    % Series is exact, X < -exp(2)
                    if X < -7.38905609893065
                        return;
                    end
                elseif r12x && Y > PI/2 && Y < 3*PI/2
                    % Region 1: near Z = -1+pi*1i, but not line of discontinuity
                    x = conj(sqrt(2*conj(Z+1-PI*1i)));                  % (20)
                    W(1) = -1+x*1i+(1440*x^2-120*x^3*1i+16*x^4 ...
                      	   +x^5*1i)/4320;                               % (21)
                elseif r12x && Y > -3*PI/2 && Y < -PI/2
                    % Region 2: near Z = -1-pi*1i, but not line of discontinuity
                    x = conj(sqrt(2*conj(Z+1+PI*1i)));                  % (22)
                  	W(1) = -1-x*1i+(1440*x^2+120*x^3*1i+16*x^4 ...
                           -x^5*1i)/4320;                               % (23)
                elseif r12x && abs(Y) <= PI/2 || ~xlteq1 && (X-1)^2+Y^2 <= PI^2
                    % Region 4: series about Z = 1
                    x = Z-1;
                    W(1) = 1+x*(1/2+x*(1/16-x*(1/192+x*(1/3072 ...
                           -(13/61440)*x))));                         	% (29)
                elseif ~xgtm2 && Y > PI && Y-PI <= -(3/4)*(X+1)
                    % Region 5: negative log series about t = Z-pi*1i
                    t = Z-pi*1i;                                       	% (26)
                    x = log(-t);
                    lzi = x/t;
                    W(1) = t-x+lzi*(1+lzi*(x/2-1+lzi*((x/3-3/2)*x+1)));	% (28)
                elseif ~xgtm2 && Y < -PI && Y+PI >= (3/4)*(X+1)
                    % Region 6: negative log series about t = Z+pi*1i
                    t = Z+PI*1i;                                        % (26)
                    x = log(-t);
                    lzi = x/t;
                    W(1) = t-x+lzi*(1+lzi*(x/2-1+lzi*((x/3-3/2)*x+1)));	% (28)
                else
                    % Region 7: log series about Z = Inf
                    x = log(Z);
                    lzi = x/Z;
                    W(1) = Z-x+lzi*(1+lzi*(x/2-1+lzi*((x/3-3/2)*x+1)));	% (25)
                end
                
                % Check for regularization: adjust Z and flip sign of W
                if X <= -0.99 && abs(Y-PI) <= 1e-2
                    Z = Z-PI*1i;                                       	% (26)
                    s = -1;
                elseif X <= -0.99 && abs(Y+PI) <= 1e-2
                    Z = Z+PI*1i;                                      	% (26)
                    s = -1;
                else
                    s = 1;
                end
            end
            
            % Residual (can be zero)
            r = Z-(W+log(s*W));                                         % (14)
            
            if abs(r) > tol
                % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
                W = W*(1+(r/(1+W))*((1+W)*(1+W+(2/3)*r)...
                    -r/2)/((1+W)*(1+W+(2/3)*r)-r));                     % (15)
                
                % Test residual
                r = Z-(W+log(s*W));                                     % (14)
                
                % Second iterative improvement via FSC method, if needed
                if abs(r) > tol
                    W = W*(1+(r/(1+W))*((1+W)*(1+W+(2/3)*r)...
                        -r/2)/((1+W)*(1+W+(2/3)*r)-r));                 % (15)
                end
            end
        end
    end
else
    isRealZ = (isreal(Z) || all(Y == 0));
    if isRealZ
        W = NaN(size(Z),dataType);
    else
        W = complex(NaN(size(Z),dataType));
    end
    
    % Special values
    W(Z == 0) = 0.5671432904097838;                             % Omega constant
    W(Z == 1) = 1;
    W(Z == 1+exp(1)) = exp(1);
    W(Z > 2^59) = Z(Z > 2^59);	% W self-saturates: X > 2^59 (abs(Y) > 2^54 too)
    
    if isRealZ
        MinZ = (Z < log(eps(realmin(dataType)))-log(2));
        W(MinZ) = 0;                                                % Z -> -Inf
    else
        MinZ = (Z(:) < log(eps(realmin(dataType)))-log(2) & Y > -PI & Y <= PI);
        W(MinZ) = 0;                                                % Z -> -Inf
        W(Z == -1+PI*1i | Z == -1-PI*1i) = -1;
        W(Z == log(1/3)-1/3+PI*1i) = -1/3;
        W(Z == log(2)-2-PI*1i) = -2;
        W(Z == (1+PI/2)*1i) = 1i;
    end
    
    i = (isnan(W(:)) & ~isnan(Z(:)));
    if any(i)
        if isRealZ
            % Region 3: series about -Inf
            j = (i & X <= -2);
            if any(j)
                x = exp(X(j));
                W(j) = x.*(1-x.*(1-x.*(36-x.*(64-125*x))/24));          % (24)
                
                % Series is exact, X < -exp(2)
                if all(~i | X < -7.38905609893065)
                    return;
                end
            end
            
            % Region 7: log series about Z = Inf
            j = (~j & X > pi+1);
            if any(j)
                t = X(j);
                x = log(t);
                lzi = x./t;
                W(j) = t-x+lzi.*(1+lzi.*(x/2-1+lzi.*((x/3-3/2).*x+1)));	% (25)
            end
            
            % Region 4: series about Z = 1
            j = (~j & X > -2);
            if any(j)
                x = X(j)-1;
                W(j) = 1+x.*(1/2+x.*(1/16-x.*(1/192+x.*(1/3072 ...
                       -(13/61440)*x))));                               % (29)
            end
            
            % No regularization
            Z = real(Z);
            s = 1;
        else
            ypi = (Y == PI);
            ympi = (Y == -PI);
            xgtm2 = (X > -2);
            xlteq1 = (X <= 1);
            r12x = (i & xgtm2 & xlteq1);
            
            % Region 3b: series about -Inf, upper line of discontinuity
            c = (i & ~xgtm2 & ypi);
            if any(c)
                x = exp(X(c));
                W(c) = ((((-125*x-64).*x-36).*x/24-1).*x-1).*x;         % (24)
                
                % Series is exact, X < -exp(2)
                i = (i & X >= -7.38905609893065);
                if all(~i)
                    return;
                end
            end
            c7 = c;
            
            % Region 3: series about -Inf, Z between lines of discontinuity
            c = (i & ~xgtm2 & Y > -PI & Y < PI);
            if any(c)
                x = exp(Z(c));
                W(c) = x.*(1-x.*(1-x.*(36+x.*(64+125*x))/24));          % (24)
                
                % Series is exact, X < -exp(2)
                if all(~i | X < -7.38905609893065)
                    return;
                end
            end
            c7 = (c7 | c);
            
            % Regions 1b and 2b: near Z = -1+/-pi*1i, lines of discontinuity
            c = (i & xgtm2 & X <= -1 & (ypi | ympi));
            if any(c)
                x = sign(Y(c)).*sqrt(-2*(X(c)+1));              	% (20, 22)
                W(c) = -1+x.*(1-x.*(1440-x.*(120+x.*(16+x)))/4320);	% (21, 23)
            end
            c7 = (c7 | c);
            
            % Region 1: near Z = -1+pi*1i, but not upper line of discontinuity
            c = (~ypi & r12x & Y > PI/2 & Y < 3*PI/2);
            if any(c)
                x = conj(sqrt(2*conj(Z(c)+1-PI*1i)));                   % (20)
                W(c) = -1+x*1i+(1440*x.^2-120*x.^3*1i+16*x.^4 ...
                       +x.^5*1i)/4320;                              	% (21)
            end
            c7 = (c7 | c);
            
            % Region 2: near Z = -1-pi*1i, but not lower line of discontinuity
            c = (~ympi & r12x & Y > -3*PI/2 & Y < -PI/2);
            if any(c)
                x = conj(sqrt(2*conj(Z(c)+1+PI*1i)));                   % (22)
                W(c) = -1-x*1i+(1440*x.^2+120*x.^3*1i+16*x.^4 ...
                       -x.^5*1i)/4320;                                	% (23)
            end
            c7 = (c7 | c);
            
            % Region 4: series about Z = 1
            c = (r12x & abs(Y) <= PI/2 | ~xlteq1 & (X-1).^2+Y.^2 <= PI^2);
            if any(c)
                x = Z(c)-1;
                W(c) = 1+x.*(1/2+x.*(1/16-x.*(1/192+x.*(1/3072 ...
                       -(13/61440)*x))));                            	% (29)
            end
            c7 = (c7 | c);
            
            % Region 5: negative log series about t = Z-pi*1i
            c = (i & ~xgtm2 & Y > PI & Y-PI <= -(3/4)*(X+1));
            if any(c)
                t = Z(c)-PI*1i;                                         % (26)
                x = log(-t);
                lzi = x./t;
                W(c) = t-x+lzi.*(1+lzi.*(x/2-1+lzi.*((x/3-3/2).*x+1)));	% (28)
            end
            c7 = (c7 | c);
            
            % Region 6b: negative log series about t = Z+pi*1i
            c = (i & ~xgtm2 & ympi);
            if any(c)
                t = X(c);                                               % (26)
                x = log(-t);
                lzi = x./t;
                W(c) = t-x+lzi.*(1+lzi.*(x/2-1+lzi.*((x/3-3/2).*x+1)));	% (28)
            end
            c7 = (c7 | c);
            
            % Region 6: negative log series about t = Z+pi*1i
            c = (i & ~xgtm2 & Y < -PI & Y+PI >= (3/4)*(X+1));
            if any(c)
                t = Z(c)+PI*1i;                                         % (26)
                x = log(-t);
                lzi = x./t;
                W(c) = t-x+lzi.*(1+lzi.*(x/2-1+lzi.*((x/3-3/2).*x+1)));	% (28)
            end
            c7 = (c7 | c);
            
            % Region 7: log series about Z = Inf
            c = (i & ~(c7 | c));
            if any(c)
                t = Z(c);
                x = log(t);
                lzi = x./t;
                W(c) = t-x+lzi.*(1+lzi.*(x/2-1+lzi.*((x/3-3/2).*x+1)));	% (25)
            end
            
            % Check for regularization: adjust Z and flip sign of W
            s1 = (X <= -0.99 & abs(Y-PI) <= 1e-2);
            Z(s1) = Z(s1)-PI*1i;                                        % (26)
            
            s2 = (X <= -0.99 & abs(Y+PI) <= 1e-2);
            Z(s2) = Z(s2)+PI*1i;                                        % (26)
            
            s = ones(size(W));
            s(s1 | s2) = -1;
        end
        
        % Residual (can be zero)
        r = Z-(W+log(s.*W));                                            % (14)
        r(MinZ) = 0;
        
        % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
        rc = (abs(r(:)) > tol);
        if any(rc)
            wr = W(rc);
            r = r(rc);
            
            W(rc) = wr.*(1+(r./(1+wr)).*((1+wr).*(1+wr+(2/3)*r)...
                    -r/2)./((1+wr).*(1+wr+(2/3)*r)-r));                 % (15)
            
            % Test residual
            r = Z-(W+log(s.*W));                                    	% (14)
            r(MinZ) = 0;
            
            rc = (abs(r(:)) > tol);
            if any(rc)
                wr = W(rc);
                r = r(rc);
                
                % Second iterative improvement via FSC method, if needed
                W(rc) = wr.*(1+(r./(1+wr)).*((1+wr).*(1+wr+(2/3)*r)...
                        -r/2)./((1+wr).*(1+wr+(2/3)*r)-r));             % (15)
            end
        end
    end
end

% Reconvert symbolic input
if isSym
    W = sym(W,'d');
end