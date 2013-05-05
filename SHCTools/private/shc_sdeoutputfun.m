function [OutputFUN,WSelect]=shc_sdeoutputfun(solver,tspan,a0,N,options)
%SHC_SDEOUTPUTFUN  Process output function arguments.
%
%   See also:
%       SHC_SDEARGUMENTS, SHC_SDEEVENTSFUN, SHC_SDERANDFUN, SHC_SDEGET,
%       SHC_SDEPLOT, FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 5-2-13
%   Revision: 1.2, 5-4-13


% Check for output function
OutputFUN = shc_sdeget(options,'OutputFUN',[],'flag');
if ~isempty(OutputFUN)
    if ~isa(OutputFUN,'function_handle')
        error('SHCTools:shc_sdeoutputfun:OutputFUNNotAFunctionHandle',...
              'OutputFUN, if specified, must be a function handle.  See %s.',...
              solver);
    end
    
    % Check for selected A components
    idxA = shc_sdeget(options,'OutputYSelect','yes','flag');
    if strcmp(idxA,'yes')
        ASelect = true;
        AAll = true;
    elseif ~strcmp(idxA,'no') && ~isempty(idxA)
        if ~isnumeric(idxA) || ~isreal(idxA) || ~all(isfinite(idxA)) ...
                || ~isvector(idxA)
            error('SHCTools:shc_sdeoutputfun:InvalidOutputASelect',...
                 ['OutputYSelect option must be a finite real numeric '...
                  'vector.  See %s.'],solver);
        end
        if any(idxA < 1) || any(idxA > N) || ~all(idxA == floor(idxA))
            error('SHCTools:shc_sdeoutputfun:InvalidIndexOutputASelect',...
                 ['OutputYSelect option must be a vector of integer indices '...
                  'no greater than the length of A0.  See %s.'],solver);
        end
        if any(diff(sort(idxA)) == 0)
            error('SHCTools:shc_sdeoutputfun:RepeatedIndexOutputASelect',...
                 ['OutputYSelect vector cannot contain repeated indices.  '...
                  'See %s.'],solver);
        end
        ASelect = true;
        AAll = false;
    else
        ASelect = false;
        AAll = false;
    end
    
    % Check for selected W components, create function handles
    idxW = shc_sdeget(options,'OutputWSelect','no','flag');
    if strcmp(idxW,'yes')
        if AAll
            OutputFUN = @(t,a,flag,w)OutputFUN(t,a,flag,w);
        elseif ASelect
            OutputFUN = @(t,a,flag,w)OutputFUN(t,a(idxA(~isempty(a),:)),flag,w);
        else
            OutputFUN = @(t,a,flag,w)OutputFUN(t,[],flag,w);
        end
        WSelect = true;
    elseif ~strcmp(idxW,'no') && ~isempty(idxW)
        if ~isnumeric(idxW) || ~isreal(idxW) || ~all(isfinite(idxW)) ...
                || ~isvector(idxW)
            error('SHCTools:shc_sdeoutputfun:InvalidOutputWSelect',...
                 ['OutputWSelect option must be a finite real numeric '...
                  'vector.  See %s.'],solver);
        end
        if any(idxW < 1) || any(idxW > N) || ~all(idxW == floor(idxW))
            error('SHCTools:shc_sdeoutputfun:InvalidIndexOutputWSelect',...
                 ['OutputWSelect option must be a vector of integer indices '...
                  'no greater than the length of A0.  See %s.'],solver);
        end
        if any(diff(sort(idxW)) == 0)
            error('SHCTools:shc_sdeoutputfun:RepeatedIndexOutputWSelect',...
                 ['OutputWSelect vector cannot contain repeated indices.  '...
                  'See %s.'],solver);
        end
        if AAll
            OutputFUN = @(t,a,flag,w)OutputFUN(t,a,flag,w(idxW(~isempty(w),:)));
        elseif ASelect
            OutputFUN = @(t,a,flag,w)OutputFUN(t,a(idxA(~isempty(a),:)),flag,...
                w(idxW(~isempty(w),:)));
        else
            OutputFUN = @(t,a,flag,w)OutputFUN(t,[],flag,...
                w(idxW(~isempty(w),:)));
        end
        WSelect = true;
    else
        if ASelect && ~AAll
            OutputFUN = @(t,a,flag,w)OutputFUN(t,a(idxA(~isempty(a),:)),flag);
        elseif ~YSelect
            OutputFUN = @(t,a,flag,w)OutputFUN(t,[],flag);
        else
            OutputFUN = @(t,a,flag,w)OutputFUN(t,a,flag);
        end
        WSelect = false;
    end
    
    % Initialize and check output function
    try
        if WSelect
            OutputFUN(tspan,a0(:),'init',zeros(N,1));
        else
            OutputFUN(tspan,a0(:),'init',[]);
        end
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                if WSelect
                    error('SHCTools:shc_sdeoutputfun:OutputFUNTooFewInputs4',...
                         ['OutputFUN must have at least four inputs.  See '...
                          'SHC_SDEPLOT.']);
                else
                    error('SHCTools:shc_sdeoutputfun:OutputFUNTooFewInputs3',...
                         ['OutputFUN must have at least three inputs.  See '...
                          'SHC_SDEPLOT.']);
                end
            case 'MATLAB:unassignedOutputs'
                error('SHCTools:shc_sdeoutputfun:OutputFUNUnassignedOutput',...
                     ['The first output of OutputFUN was not assigned.  See '...
                      'SHC_SDEPLOT.']);
            case 'MATLAB:minrhs'
                error('SHCTools:shc_sdeoutputfun:OutputFUNTooManyInputs',...
                     ['OutputFUN requires one or more input arguments '...
                      '(parameters) that were not supplied.  See '...
                      'SHC_SDEPLOT.']);
            otherwise
                rethrow(err);
        end
    end
else
    WSelect = false;
end