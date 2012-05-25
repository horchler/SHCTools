function y=shc_hypergeom2F1q(a,b,c,z,tol)
%SHC_HYPERGEOM2F1Q  Gaussian or ordinary hypergeometric function for ABS(Z) < 1
%
%   Note:
%       Unless, C = A or C = B, if C is an integer <=0, NaN is returned. 

%   Andrew D. Horchler, adh9@case.edu, Created 5-12-13
%   Revision: 1.0, 5-23-12


% Calculate Gaussian hypergeometric function via series summation
if isempty(z)
	y = z;
elseif abs(z) >= 1
    if abs(z) == 1
        error('SHTools:shc_hypergeom2F1q:ZSingularity',...
             ['Singularity: The standard Gaussian hypergeometric function '...
              'is only defined on ABS(Z) < 1']);
    else
        error('SHTools:shc_hypergeom2F1q:ZUndefined',...
             ['The standard Gaussian hypergeometric function is only '...
              'defined on ABS(Z) < 1']);
    end
else
    % Set relative tolerance and maximum number of iterations
    if nargin < 5
        tol = eps(superiorfloat(a,b,c,z));
    end
    itermax = 2^15;
    
    yp = 0;
    i = 0;
    if a == 1
        if b == c
            yi = z;
            y = 1+yi;
            while abs(y-yp) > tol*y && i < itermax
                yp = y;
                i = i+1;
                yi = yi*z;
                y = y+yi;
            end
        else
            yi = b*z/c;
            y = 1+yi;
            while abs(y-yp) > tol*y && i < itermax
                yp = y;
                i = i+1;
                yi = yi*(b+i)*z/(c+i);
                y = y+yi;
            end
        end
    elseif b == 1
        if a == c
            yi = z;
            y = 1+yi;
            while abs(y-yp) > tol*y && i < itermax
                yp = y;
                i = i+1;
                yi = yi*z;
                y = y+yi;
            end
        else
            yi = a*z/c;
            y = 1+yi;
            while abs(y-yp) > tol*y && i < itermax
                yp = y;
                i = i+1;
                yi = yi*(a+i)*z/(c+i);
                y = y+yi;
            end
        end
    elseif a == c
        yi = b*z;
        y = 1+yi;
        while abs(y-yp) > tol*y && i < itermax
            yp = y;
            i = i+1;
            yi = yi*(b+i)*z/(i+1);
            y = y+yi;
        end
    elseif b == c
        yi = a*z;
        y = 1+yi;
        while abs(y-yp) > tol*y && i < itermax
            yp = y;
            i = i+1;
            yi = yi*(a+i)*z/(i+1);
            y = y+yi;
        end
    else 
        yi = a*b*z/c;
        y = 1+yi;
        while abs(y-yp) > tol*y && i < itermax
            yp = y;
            i = i+1;
            yi = yi*(a+i)*(b+i)*z/((i+1)*(c+i));
            y = y+yi;
        end
    end
end