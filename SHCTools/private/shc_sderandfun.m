function [RandFUN,ResetStream]=shc_sderandfun(solver,dataType,options)
%SHC_SDERANDFUN  Process random function arguments.
%
%   See also:
%       SHC_SDEARGUMENTS, SHC_SDEEVENTSFUN, SHC_SDEOUTPUTFUN, SHC_SDEGET,
%       FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 5-2-13
%   Revision: 1.3, 11-28-13


ResetStream = [];
RandFun = shc_sdeget(options,'RandFun',[],'flag');
if isempty(RandFun) || isa(RandFun,'RandStream')
    if isempty(RandFun)
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
        Stream = RandFun;
    end
    
    % Create function handle to be used for generating Wiener increments
    RandFUN = @(M,N)randn(Stream,M,N,dataType);
elseif isa(RandFun,'function_handle')
    % Create function handle to be used for generating Wiener increments
    RandFUN = @(M,N)RandFunCheck(RandFun,M,N,solver);
elseif sde_ismatrix(RandFun)
    RandFUN = [];
else
    error('SHCTools:shc_sderandfun:InvalidRandFun',...
     	 ['RandFun must be a RandStream object, function handle, or matrix '...
          'of integrated Wiener increments.  See %s.'],solver);
end


function R=RandFunCheck(fun,M,N,solver)
try
    R = feval(fun,M,N);
    if ~shc_ismatrix(R) || isempty(R) || ~isfloat(R)
        error('SHCTools:shc_sderandfun:RandFunCheck:RandFUNNot2DArray3',...
             ['RandFUN must return a non-empty matrix of '...
              'floating-point values.  See %s.'],solver);
    end
    [m,n] = size(R);
    if m ~= M || n ~= N
        error('SHCTools:shc_sderandfun:RandFunCheck:RandFUNDimensionMismatch3',...
             ['The specified alternative RandFUN did not '...
              'output a %d by %d matrix as requested.  '...
              'See %s.'],M,N,solver);
    end
catch err
    switch err.identifier
        case 'MATLAB:TooManyInputs'
            error('SHCTools:shc_sderandfun:RandFunCheck:RandFUNTooFewInputs',...
                 ['RandFUN must have at least two inputs.  '...
                  'See %s.'],solver);
        case 'MATLAB:TooManyOutputs'
            error('SHCTools:shc_sderandfun:RandFunCheck:RandFUNNoOutput',...
                 ['The output of RandFUN was not specified. '...
                  'RandFUN must return a non-empty matrix.  '...
                  'See %s.'],solver);
        case 'MATLAB:unassignedOutputs'
            error('SHCTools:shc_sderandfun:RandFunCheck:RandFUNUnassignedOutput',...
                 ['The first output of RandFUN was not '...
                  'assigned.  See %s.'],solver);
        case 'MATLAB:minrhs'
            error('SHCTools:shc_sderandfun:RandFunCheck:RandFUNTooManyInputs',...
                 ['RandFUN must not require more than two '...
                  'inputs.  See %s.'],solver);
        otherwise
            rethrow(err);
    end
end

%{

% Check if alternative random number generator function or W matrix specified
isNoRandFUN = (~isfield(options,'RandFUN') || isempty(options.RandFUN));
if isNoRandFUN || isa(options.RandFUN,'RandStream')
    if isNoRandFUN
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
        end
    else
        % User-specified RandStream object, ignore RandSeed and Antithetic
        Stream = options.RandFUN;
    end
    
    % Create function handle to be used for generating Wiener increments
    RandFUN = @(M,N)randn(Stream,M,N,dataType);
    
    % Function to be call on completion or early termination of integration
    ResetStream = onCleanup(@()shc_sdereset_stream(Stream));
else
    % Function handle for generating Wiener increments created in main function
    RandFUN = [];
    ResetStream = [];
end
%}