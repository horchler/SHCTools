function shc_sdereset_stream(Stream)
%SHC_SDERESET_STREAM  Reset global stream or delete custom random stream.
%
%   See also:
%       SHC_LV_INTEGRATE, SHC_SDERANDFUN, ONCLEANUP, FUNCTION_HANDLE, RANDSTREAM
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 12-30-12
%   Revision: 1.2, 5-4-13


% Called by onCleanup to reset antihetic property for global stream
try
    isGlobal = isequal(Stream,RandStream.getGlobalStream);
catch                                                       %#ok<CTCH>
    isGlobal = isequal(Stream,RandStream.getDefaultStream);	%#ok<GETRS>
end
if isGlobal
    set(Stream,'Antithetic',~Stream.Antithetic);
else
    delete(Stream);
end