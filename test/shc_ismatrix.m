function bool=shc_ismatrix(V)
%SHC_ISMATRIX  Replicate functionality of builtin ismatrix for pre-R2010b Matlab
%
%

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-9-12
%   Revision: 1.0, 6-30-12


try
    bool = ismatrix(V);
catch err
    if strcmp(err.identifier,'MATLAB:UndefinedFunction')
        bool = ndims(V) == 2 || size(V(:,:,:),3) == 1;	%#ok<ISMAT>
    else
        rethrow(err);
    end
end