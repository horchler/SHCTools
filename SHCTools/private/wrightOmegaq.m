function W=wrightOmegaq(Z)
%WRIGHTOMEGAQ  Wright omega function, a solution of the equation Y+LOG(Y) = Z.
%   W = WRIGHTOMEGAQ(Z) performs double-precision evaluation of the Wright omega
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
%       absolute sense) to Z, e.g., Z <= -714.84989998498+B*1i, -pi < B <= pi.
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

% 	The inverse Wright omega function is defined as (Corless & Jeffrey, 2002):
%       Y+LOG(Y)+2*pi*1i,   -Inf < Y < -1
%       -1 +/- pi*1i,       Y = 1
%       Y+LOG(Y),           otherwise

%   WRIGHTOMEGAQ is 3 to 4 orders of magnitude faster than WRIGHTOMEGA for
%   double-precision arrays. Additionally, it has less numeric error, properly
%   evaluates values greater than 2^28, supports single-precision evaluation,
%   and handles NaN inputs.

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-12-12
%   Revision: 1.0, 2-11-12


% Convert symbolic input, converted back at end
isSym = isa(Z,'sym') || ischar(Z);
if isSym
    try
        Z = double(sym(Z));
    catch ME
        if strcmp(ME.identifier,'symbolic:sym:double:cantconvert')
            error('SHCTools:wrightOmegaq:InvalidSymbolicZ',...
                 ['Symbolic array must contain only numeric values and no '...
                  'expressions containing variables.'])
        else
            rethrow(ME);
        end
    end
elseif ~isfloat(Z)
    error('SHCTools:wrightOmegaq:InvalidZ',...
         ['Input Z must be a floating point array or an array of symbolic '...
          'values.']);
end

dataType = class(Z);
tol = eps(dataType);

% Special scalar-optimized case, W(1) used in order retain datatype
if isscalar(Z)
    % Special values
    if isempty(Z) || isnan(Z) || Z > 2^59
        W = Z;  % W self-saturates if real(Z) > 2^59 (abs(imag(Z)) > 2^54 also)
    elseif Z == 0
        W(1) = 0.5671432904097838;                 	% Omega constant
    elseif Z == 1
        W(1) = 1;
    elseif Z == 1+exp(1)
        W(1) = exp(1);
    elseif isreal(Z) || imag(Z == 0)
        if Z < log(eps(realmin(dataType)))-log(2)
            W(1) = 0;                             	% Z -> -Inf
        else
            % Real case
            Z = real(Z);

            if Z <= -2
                % Region 3: series about -Inf
                ez = exp(Z);
                w = ez*(1-ez*(1-ez*(3/2+ez*(8/3+(125/24)*ez))));        % (24)
            elseif Z > pi+1
                % Region 7: log series about z = Inf
                lz = log(Z);
                lzi = lz/Z;
                w = Z-lz+lzi*(1+lzi*(lz/2-1+lzi*((lz/3-3/2)*lz+1)));    % (25)
            else
                % Region 4: series about z = 1
                zc = Z-1;
                w = 1+zc*(1/2+zc*(1/16-zc*(1/192+zc*(1/3072 ...
                    -(13/61440)*zc))));                                 % (29)
            end
            
            % Residual
            r = Z-(w+log(w));                                           % (14)
            
            % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
            w = w*(1+(r/(1+w))*((1+w)*(1+w+(2/3)*r)...
                -r/2)/((1+w)*(1+w+(2/3)*r)-r));                         % (15)
            
            % Test residual
            r = Z-(w+log(w));                                           % (14)
            
            % Second iterative improvement via FSC method, if needed
            if abs(r) > tol
                w = w*(1+(r/(1+w))*((1+w)*(1+w+(2/3)*r)...
                    -r/2)/((1+w)*(1+w+(2/3)*r)-r));                     % (15)
            end
            
            W(1) = w;
        end
    else
        % Special complex values
        if Z < log(eps(realmin(dataType)))-log(2) && imag(Z) > -pi ...
                && imag(Z) <= pi
            W(1) = 0;                               	% Z -> -Inf
        elseif Z == -1+pi*1i || Z == -1-pi*1i
            W(1) = -1;
        elseif Z == log(1/3)-1/3+pi*1i
            W(1) = -1/3;
        elseif Z == log(2)-2-pi*1i
            W(1) = -2;
        elseif Z == (1+pi/2)*1i
            W(1) = 1i;
        else
            % Complex Case
            x = real(Z);
            y = imag(Z);
            
            xgtm2 = (x > -2);
            xlteq1 = (x <= 1);
            r12x = (xgtm2 && xlteq1);
            
            if r12x && y > pi/2 && y < 3*pi/2
                % Region 1: near z = -1+pi*1i
                zc = conj(sqrt(2*conj(Z+1-pi*1i)));                     % (20)
                w = -1+zc*1i+zc^2/3-zc^3*1i/36+zc^4/270 ...
                    +zc^5*1i/4320;                                      % (21)
            elseif r12x && y > -3*pi/2 && y < -pi/2
                % Region 2: near z = -1-pi*1i
                zc = conj(sqrt(2*conj(Z+1+pi*1i)));                     % (22)
                w = -1-zc*1i+zc^2/3+zc^3*1i/36+zc^4/270 ...
                    -zc^5*1i/4320;                                      % (23)
            elseif ~xgtm2 && y > -pi && y <= pi
                % Region 3: series about -Inf, z between lines of discontinuity
                if y == pi
                    ez = -exp(x);	% Avoid small imaginary error, common case 
                else
                    ez = exp(Z);
                end
                w = ez*(1-ez*(1-ez*(3/2+ez*(8/3+(125/24)*ez))));        % (24)
            elseif xgtm2 && xlteq1 && abs(y) <= pi/2 || ~xlteq1 ...
                    && (x-1)^2+y^2 <= pi^2
                % Region 4: series about z = 1
                z1 = Z-1;
                w = 1+z1*(1/2+z1*(1/16-z1*(1/192+z1*(1/3072 ...
                    -(13/61440)*z1))));                                 % (29)
            elseif ~xgtm2 && y > pi && y-pi <= -(3/4)*(x+1)
                % Region 5: negative log series about t = z-pi*1i
                t = Z-pi*1i;                                            % (26)
                lz = log(-t);
                lzi = lz/t;
                w = t-lz+lzi*(1+lzi*(lz/2-1+lzi*((lz/3-3/2)*lz+1)));    % (28)
            elseif ~xgtm2 && y <= -pi && y+pi >= (3/4)*(x+1)
                % Region 6: negative log series about t = z+pi*1i
                t = Z+pi*1i;                                            % (26)
                lz = log(-t);
                lzi = lz/t;
                w = t-lz+lzi*(1+lzi*(lz/2-1+lzi*((lz/3-3/2)*lz+1)));    % (28)
            else
                % Region 7: log series about z = Inf
                lz = log(Z);
                lzi = lz/Z;
                w = Z-lz+lzi*(1+lzi*(lz/2-1+lzi*((lz/3-3/2)*lz+1)));    % (25)
            end
            
            % Check for regularization: adjust z and flip sign of w
            if x <= -0.99 && abs(y-pi) <= 1e-2
                Z = Z-pi*1i;                                            % (26)
                s = -1;
            elseif x <= -0.99 && abs(y+pi) <= 1e-2
                Z = Z+pi*1i;                                            % (26)
                s = -1;
            else
                s = 1;
            end
            
            % Residual (can be zero)
            r = Z-(w+log(s*w));                                         % (14)
            
            if abs(r) > tol
                % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
                w = w*(1+(r/(1+w))*((1+w)*(1+w+(2/3)*r)...
                    -r/2)/((1+w)*(1+w+(2/3)*r)-r));                     % (15)
                
                % Test residual
                r = Z-(w+log(s*w));                                     % (14)
                
                % Second iterative improvement via FSC method, if needed
                if abs(r) > tol
                    w = w*(1+(r/(1+w))*((1+w)*(1+w+(2/3)*r)...
                        -r/2)/((1+w)*(1+w+(2/3)*r)-r));                 % (15)
                end
            end
            
            W(1) = w;
        end
    end
else
    % General array input vectorized case
    isRealZ = (isreal(Z) || all(imag(Z(:)) == 0));
    if isRealZ
        W = NaN(size(Z),dataType);
    else
        W = complex(NaN(size(Z),dataType));
    end
    if ~isempty(W)
        % Special values
        W(Z == 0) = 0.5671432904097838;                     % Omega constant
        W(Z == 1) = 1;
        W(Z == 1+exp(1)) = exp(1);
        
        % W self-saturates if real(Z) > 2^59 (abs(imag(Z)) > 2^54 also)
        W(Z > 2^59) = Z(Z > 2^59);
        if isRealZ
            W(Z < log(eps(realmin(dataType)))-log(2)) = 0;	% Z -> -Inf
        else
            W(Z < log(eps(realmin(dataType)))-log(2) & imag(Z) > -pi ...
                & imag(Z) <= pi) = 0;                       % Z -> -Inf
            W(Z == -1+pi*1i | Z == -1-pi*1i) = -1;
            W(Z == log(1/3)-1/3+pi*1i) = -1/3;
            W(Z == log(2)-2-pi*1i) = -2;
            W(Z == (1+pi/2)*1i) = 1i;
        end
        
        iz = (isnan(W(:)) & ~isnan(Z(:)));
        if any(iz)
            z = Z(iz);
            w = W(iz);
            
            % Numbers in parentheses refer to equations in Lawrence, et al. 2012
            if isRealZ
                % Real case
                z = real(z);
                
                % Region 3: series about -Inf
                c = (z <= -2);
                if any(c)
                    zc = exp(z(c));
                    w(c) = zc.*(1-zc.*(1-zc.*(3/2+zc.*(8/3 ...
                        +(125/24)*zc))));                               % (24)
                end
                c4 = c;
                
                % Region 7: log series about z = Inf
                c = (z > pi+1);
                if any(c)
                    zc = z(c);
                    lz = log(zc);
                    lzi = lz./zc;
                    w(c) = zc-lz+lzi.*(1+lzi.*(lz/2-1 ...
                        +lzi.*((lz/3-3/2).*lz+1)));                     % (25)
                end
                
                % Region 4: series about z = 1
                c = ~(c4 | c);
                if any(c)
                    zc = z(c)-1;
                    w(c) = 1+zc.*(1/2+zc.*(1/16-zc.*(1/192+zc.*(1/3072 ...
                        -(13/61440)*zc))));                           	% (29)
                end
                
                % Residual
                r = z-(w+log(w));                                     	% (14)
                
                % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
                w = w.*(1+(r./(1+w)).*((1+w).*(1+w+(2/3)*r)...
                    -r/2)./((1+w).*(1+w+(2/3)*r)-r));                  	% (15)
                
                % Test residual
                r = z-(w+log(w));                                      	% (14)
                
                rc = (abs(r) > tol);
                if any(rc)
                    wr = w(rc);
                    r = r(rc);
                    
                    % Second iterative improvement via FSC method, if needed
                    w(rc) = wr.*(1+(r./(1+wr)).*((1+wr).*(1+wr+(2/3)*r)...
                        -r/2)./((1+wr).*(1+wr+(2/3)*r)-r));            	% (15)
                end
            else
                % Complex Case
                x = real(z);
                y = imag(z);
                
                xgtm2 = (x > -2);
                xlteq1 = x <= 1;
                r12x = (xgtm2 & xlteq1);
                
                % Region 1: near z = -1+pi*1i
                c = (r12x & y > pi/2 & y < 3*pi/2);
                if any(c)
                    zc = conj(sqrt(2*conj(z(c)+1-pi*1i)));            	% (20)
                    w(c) = -1+zc*1i+zc.^2/3-zc.^3*1i/36+zc.^4/270 ...
                        +zc.^5*1i/4320;                                	% (21)
                end
                c7 = c;
                
                % Region 2: near z = -1-pi*1i
                c = (r12x & y > -3*pi/2 & y < -pi/2);
                if any(c)
                    zc = conj(sqrt(2*conj(z(c)+1+pi*1i)));            	% (22)
                    w(c) = -1-zc*1i+zc.^2/3+zc.^3*1i/36+zc.^4/270 ...
                        -zc.^5*1i/4320;                                	% (23)
                end
                c7 = (c7 | c);
                
                % Region 3: series about -Inf, z between lines of discontinuity
                c = (~xgtm2 & y > -pi & y <= pi);
                if any(c)
                    ypi = (y == pi);
                    if any(ypi)     % Avoid small imaginary error, common case
                        zc = exp(y(c));
                        zc(c & ypi) = -1;
                        zc = zc.*exp(x(c));
                    else
                        zc = exp(z(c));
                    end
                    w(c) = zc.*(1-zc.*(1-zc.*(3/2+zc.*(8/3 ...
                        +(125/24)*zc))));                               % (24)
                end
                c7 = (c7 | c);
                
                % Region 4: series about z = 1
                c = (xgtm2 & xlteq1 & abs(y) <= pi/2 | ~xlteq1 ...
                    & (x-1).^2+y.^2 <= pi^2);
                if any(c)
                    zc = z(c)-1;
                    w(c) = 1+zc.*(1/2+zc.*(1/16-zc.*(1/192+zc.*(1/3072 ...
                        -(13/61440)*zc))));                            	% (29)
                end
                c7 = (c7 | c);
                
                % Region 5: negative log series about t = z-pi*1i
                c = (~xgtm2 & y > pi & y-pi <= -(3/4)*(x+1));
                if any(c)
                    zc = z(c)-pi*1i;                                  	% (26)
                    lz = log(-zc);
                    lzi = lz./zc;
                    w(c) = zc-lz+lzi.*(1+lzi.*(lz/2-1 ...
                        +lzi.*((lz/3-3/2).*lz+1)));                   	% (28)
                end
                c7 = (c7 | c);
                
                % Region 6: negative log series about t = z+pi*1i
                c = (~xgtm2 & y <= -pi & y+pi >= (3/4)*(x+1));
                if any(c)
                    zc = z(c)+pi*1i;                                 	% (26)
                    lz = log(-zc);
                    lzi = lz./zc;
                    w(c) = zc-lz+lzi.*(1+lzi.*(lz/2-1 ...
                        +lzi.*((lz/3-3/2).*lz+1)));                  	% (28)
                end
                
                % Region 7: log series about z = Inf
                c = ~(c7 | c);
                if any(c)
                    zc = z(c);
                    lz = log(zc);
                    lzi = lz./zc;
                    w(c) = zc-lz+lzi.*(1+lzi.*(lz/2-1 ...
                        +lzi.*((lz/3-3/2).*lz+1)));                   	% (25)
                end
                
                % Check for regularization: adjust z and flip sign of w
                s1 = (x <= -0.99 & abs(y-pi) <= 1e-2);
                z(s1) = z(s1)-pi*1i;                                  	% (26)
                
                s2 = (x <= -0.99 & abs(y+pi) <= 1e-2);
                z(s2) = z(s2)+pi*1i;                                 	% (26)
                
                s = double(~(s1 | s2));
                s(s == 0) = -1;
                
                % Residual (can be zero)
                r = z-(w+log(s.*w));                                  	% (14)
                
                % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
                rc = (abs(r) > tol);
                if any(rc)
                    wr = w(rc);
                    r = r(rc);
                    
                    w(rc) = wr.*(1+(r./(1+wr)).*((1+wr).*(1+wr+(2/3)*r)...
                        -r/2)./((1+wr).*(1+wr+(2/3)*r)-r));            	% (15)
                    
                    % Test residual
                    r = z-(w+log(s.*w));                               	% (14)
                    
                    rc = (abs(r) > tol);
                    if any(rc)
                        wr = w(rc);
                        r = r(rc);
                        
                        % Second iterative improvement via FSC method, if needed
                        w(rc) = wr.*(1+(r./(1+wr)).*((1+wr).*(1+wr+(2/3)*r)...
                            -r/2)./((1+wr).*(1+wr+(2/3)*r)-r));        	% (15)
                    end
                end
            end

            W(iz) = w;
        end
    end
end

% Set any all zero imaginary outputs to real
if all(imag(W) == 0)
    W = real(W);
end

% Reconvert symbolic input
if isSym
    W = sym(W,'d');
end