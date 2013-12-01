classdef CACHE < handle
    %CACHE  
    %
    %   See also CACHE/IN, CACHE/OUT, CACHE/FIND, CACHE/COMPARE, CACHE/RESIZE,
    %       CACHE/COMPARISON, CACHE/CLEAR
    
    %   Andrew D. Horchler, adh9 @ case . edu, Created 11-29-13
    %   Revision: 1.0, 12-1-13
    
    
    properties (SetAccess=private)
        Size
        Used = 0;
    end
    
    properties (SetAccess=private,Hidden)
        numin
        numout
        inputs
        outputs
        tests = {{}};
    end
    
    % --------------------------------------------------------------------------
    
    methods
        
        function obj = CACHE(Size,varargin)
            if nargin == 0
                Size = 1;
            elseif ~CACHE.isIndex(Size) || Size < 1
                error('CACHE:InvalidCacheSize',...
                     ['Cache size must be a finite real integer greater '...
                      'than zero.']);
            end
            obj.Size = Size;
            obj.numin = zeros(Size,1);
            obj.numout = zeros(Size,1);
            obj.inputs = cell(Size,1);
            obj.outputs = cell(Size,1);
            if nargin > 1
                IN(obj,[],varargin{:});
            end
        end
        
        function idx = IN(obj,idx,varargin)
            if nargin > 1
                if isempty(idx)
                    idx = FIND(obj,varargin{:});
                    if isempty(idx)
                        if obj.Used < obj.Size
                            obj.Used = obj.Used+1;
                            obj.numin(obj.Used) = nargin-2;
                            obj.inputs{obj.Used} = varargin;
                        else
                            obj.numin = [obj.numin(2:end);nargin-2];
                            [obj.inputs{:}] = deal(obj.inputs{2:end},varargin);
                        end
                    end
                else
                    if ~CACHE.isIndex(idx) || idx < 1 || idx > obj.Used
                        error('CACHE:IN:InvalidCacheIndex',...
                             ['Cache index must be a finite real integer '...
                              'greater than zero and less than or equal to '...
                              'the used cache.']);
                    end
                    obj.numin(idx) = nargin-2;
                    obj.inputs{idx} = varargin;
                end
            else
                idx = obj.Used;
            end
        end
        
        function [n,varargout] = OUT(obj,idx,varargin)
            if obj.Used == 0
                error('CACHE:OUT:NotInitialized',...
                   	  'Cache has not been initialized.');
            elseif nargin == 1 || isempty(idx)
             	idx = obj.Used;
        	elseif ~CACHE.isIndex(idx) || idx < 1 || idx > obj.Used
                error('CACHE:OUT:InvalidCacheIndex',...
                     ['Cache index must be a finite real integer greater '...
                      'than zero and less than or equal to the used cache.']);
            end
            
            if nargin > 2
                obj.numout(idx) = nargin-2;
             	obj.outputs{idx} = varargin;
            end
            
            n = obj.numout(idx);
            if nargout > 1
                if nargout-1 > n
                    if n == 0
                        error('CACHE:OUT:NoOutputs',...
                              'No outputs have been set for cache index %d.',...
                              idx);
                    else
                        error('CACHE:OUT:TooManyOutputs',...
                              'Too many output arguments.');
                    end
                end
                v = 1:min(nargout-1,n);
                [varargout{v}] = obj.outputs{idx}{v};
            end
        end
        
        function idx = FIND(obj,varargin)
            idx = [];
            if ~isempty(varargin)
                for i = 1:obj.Used
                    if isscalar(obj.tests)
                        testcell = obj.tests{1};
                    else
                        testcell = obj.tests{i};
                    end
                    if isempty(testcell)
                        if nargin-1 == obj.numin(i) ...
                                && isequaln(obj.inputs{i},varargin)
                            idx = i;
                            break;
                        end
                    else
                        tf = (nargin-1 == obj.numin(i));
                        j = 1;
                        while tf && j <= obj.numin(i)
                            if isequal(size(obj.inputs{i}{j}),size(varargin{j}))
                                if isscalar(testcell)
                                    test = testcell{1};
                                else
                                    test = testcell{j};
                                end
                            
                                tf = all(feval(test,varargin{j},...
                                                    obj.inputs{i}{j}));
                                j = j+1;
                            else
                                tf = false;
                            end
                        end
                        if tf
                            idx = i;
                            break;
                        end
                    end
                end
            end
        end
        
        function tf = COMPARE(obj,varargin)
            tf = ~isempty(FIND(obj,varargin{:}));
        end
        
        function obj = COMPARISON(obj,tests,idx)
            if nargin > 1
                if nargin > 2
                    if isempty(idx)
                        idx = obj.Used;
                    elseif ~CACHE.isIndex(idx) || idx < 1 || idx > obj.Size
                        error('CACHE:COMPARISON:InvalidCacheIndex',...
                             ['Cache index must be a finite real integer '...
                              'greater than zero and less than or equal to '...
                              'the cache size.']);
                    end
                else
                    idx = 1:obj.Used;
                end
            
                if iscell(tests)
                    if ~(isscalar(tests) || isvector(tests) ...
                            && all(length(tests) == obj.numin(idx)))
                        error('CACHE:COMPARISON:DimensionMismatch',...
                             ['Cell arrays of comparison test functions '...
                              'must either be scalars or the same length as '...
                              'the number of inputs for the corresponding '...
                              'cache index.']);
                    end
                    if ~all(cellfun('isclass',tests,'function_handle')) ...
                            || any(cellfun(@(x)any(nargin(x) == [0 1]),...
                            tests)) || any(cellfun(@nargout,tests) == 0)
                        error('CACHE:COMPARISON:InvalidTestCell',...
                             ['Comparison test functions must be function '...
                              'handles that accept at least two inputs and '...
                              'return at least one ouput.']);
                    end
                else
                    if ~isscalar(tests) || ~isa(tests,'function_handle') ...
                            || any(nargin(tests) == [0 1]) ...
                            || nargout(tests) == 0
                        error('CACHE:COMPARISON:InvalidTest',...
                             ['Comparison test functions must be function '...
                              'handles that accept at least two inputs and '...
                              'return at least one ouput.']);
                    end
                    tests = {tests};
                end
                
                if nargin > 2
                    if isscalar(obj.tests)
                        obj.tests = repmat(obj.tests,obj.Size,1);
                    end
                    obj.tests{idx} = tests;
                else
                    obj.tests = {tests};
                end
            end
        end
        
        function obj = RESIZE(obj,Size)
            if nargin == 1
                Size = 1;
            elseif ~CACHE.isIndex(Size) || Size < 0
                error('CACHE:RESIZE:InvalidCacheSize',...
                     ['Resized cache size must be a finite real integer '...
                      'greater than or equal to zero.']);
            end
            if Size > obj.Size
                obj.Size = Size;
                obj.numin = [obj.numin;zeros(Size-obj.Size,1)];
                obj.numout = [obj.numout;zeros(Size-obj.Size,1)];
            	obj.inputs = [obj.inputs;cell(Size-obj.Size,1)];
                obj.outputs = [obj.outputs;cell(Size-obj.Size,1)];
                if ~isscalar(obj.tests)
                    obj.tests = [obj.tests;repmat({{}},Size-obj.Size,1)];
                end
            elseif Size < obj.Size
                if Size == 0
                    obj = CLEAR(obj);
                else
                    obj.Size = Size;
                    obj.Used = min(obj.Used,Size);
                    obj.numin = obj.numin(end-Size+1:end);
                    obj.numout = obj.numout(end-Size+1:end);
                    obj.inputs = obj.inputs(end-Size+1:end);
                    obj.outputs = obj.outputs(end-Size+1:end);
                    if ~isscalar(obj.tests)
                        obj.tests = obj.tests(end-Size+1:end);
                    end
                end
            end
        end
        
        function obj = CLEAR(obj)
            obj = delete(obj);
        end
        
        function obj = delete(obj)
            obj.Size = 0;
            obj.Used = 0;
            obj.numin = [];
            obj.numout = [];
          	obj.inputs = [];
            obj.outputs = [];
            obj.tests = {{}};
        end
        
    end
    
    % --------------------------------------------------------------------------
    
    methods (Access=private,Static,Hidden)
        
        function tf = isIndex(idx)
            tf = (isscalar(idx) && isnumeric(idx) && isreal(idx) ...
                && isfinite(idx) && idx == floor(idx));
        end
        
    end
    
end