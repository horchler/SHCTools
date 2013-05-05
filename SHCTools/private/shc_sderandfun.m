function [RandFUN,ResetStream]=shc_sderandfun(solver,dataType,options)
%SHC_SDERANDFUN  Process random function arguments.
%
%   See also:
%       SHC_SDEARGUMENTS, SHC_SDEEVENTSFUN, SHC_SDEOUTPUTFUN, SHC_SDEGET,
%       FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 5-2-13
%   Revision: 1.2, 5-4-13


% Check if alternative random number generator function or W matrix specified
isNoRandFUN = (isempty(options) || ~isfield(options,'RandFUN'));
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