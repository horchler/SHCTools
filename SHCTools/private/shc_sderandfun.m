function [RandFUN,ResetStream]=shc_sderandfun(solver,dataType,options)
%SHC_SDERANDFUN  Process random function arguments.
%
%   See also:
%       SHC_SDEARGUMENTS, SHC_SDEEVENTSFUN, SHC_SDEOUTPUTFUN, SHC_SDEGET,
%       FUNCTION_HANDLE
        
%   Andrew D. Horchler, horchler @ gmail . com, Created 5-2-13
%   Revision: 1.3, 11-28-13


ResetStream = [];
RandFUN = shc_sdeget(options,'RandFUN',[],'flag');
if isempty(RandFUN) || isa(RandFUN,'RandStream')
    if isempty(RandFUN)
        % Use Matlab's random number generator for normal variates
        RandSeed = shc_sdeget(options,'RandSeed',[],'flag');
        if ~isempty(RandSeed)
            if ~isscalar(RandSeed) || ~isnumeric(RandSeed) ...
                    || ~isreal(RandSeed) || ~isfinite(RandSeed) ...
                    || RandSeed >= 2^32 || RandSeed < 0
                error('SHCTools:shc_sderandfun:InvalidRandSeed',...
                     ['RandSeed must be a non-negative integer value less '...
                      'than 2^32.  See %s.'],solver);
            end
            % Create new stream based on seed value
            Stream = RandStream.create('mt19937ar','Seed',RandSeed);
        else
            % Use default stream
            try
                Stream = RandStream.getGlobalStream;
            catch                                       %#ok<CTCH>
                Stream = RandStream.getDefaultStream;	%#ok<GETRS>
            end
        end

        % Set property if antithetic random variates option is specified
        Antithetic = strcmp(shc_sdeget(options,'Antithetic','no','flag'),'yes');
        if Antithetic ~= Stream.Antithetic
            set(Stream,'Antithetic',Antithetic);
            
            % Reset property on completion or early termination of integration
            ResetStream = onCleanup(@()set(Stream,'Antithetic',...
                                           ~Stream.Antithetic));
        end
    else
        % User-specified RandStream object, ignore RandSeed and Antithetic
        Stream = RandFUN;
    end
    
    % Create function handle to be used for generating Wiener increments
    RandFUN = @(M,N)randn(Stream,M,N,dataType);
elseif isa(RandFUN,'function_handle')
    % Create function handle to be used for generating Wiener increments
    RandFUN = @(M,N)RandFUNCheck(RandFUN,M,N,solver);
elseif shc_ismatrix(RandFUN)
    RandFUN = [];
else
    error('SHCTools:shc_sderandfun:InvalidRandFUN',...
     	 ['RandFUN must be a RandStream object, function handle, or matrix '...
          'of integrated Wiener increments.  See %s.'],solver);
end


function R=RandFUNCheck(fun,M,N,solver)
try
    R = feval(fun,M,N);
    if ~shc_ismatrix(R) || isempty(R) || ~isfloat(R)
        error('SHCTools:shc_sderandfun:RandFUNCheck:RandFUNNot2DArray3',...
             ['RandFUN must return a non-empty matrix of floating-point '...
              'values.  See %s.'],solver);
    end
    [m,n] = size(R);
    if m ~= M || n ~= N
        error('SHCTools:shc_sderandfun:RandFUNCheck:RandFUNDimensionMismatch3',...
             ['The specified alternative RandFUN did not output a %d-by-%d '...
              'matrix as requested.  See %s.'],M,N,solver);
    end
catch err
    switch err.identifier
        case 'MATLAB:TooManyInputs'
            error('SHCTools:shc_sderandfun:RandFUNCheck:RandFUNTooFewInputs',...
                  'RandFUN must have at least two inputs.  See %s.',solver);
        case 'MATLAB:TooManyOutputs'
            error('SHCTools:shc_sderandfun:RandFUNCheck:RandFUNNoOutput',...
                 ['The output of RandFUN was not specified. RandFUN must '...
                  'return a non-empty matrix.  See %s.'],solver);
        case 'MATLAB:unassignedOutputs'
            error('SHCTools:shc_sderandfun:RandFUNCheck:RandFUNUnassignedOutput',...
                 ['The first output of RandFUN was not '...
                  'assigned.  See %s.'],solver);
        case 'MATLAB:minrhs'
            error('SHCTools:shc_sderandfun:RandFUNCheck:RandFUNTooManyInputs',...
                 ['RandFUN must not require more than two '...
                  'inputs.  See %s.'],solver);
        otherwise
            rethrow(err);
    end
end