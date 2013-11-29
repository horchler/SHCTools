classdef CACHE < handle
    %CACHE  
    %
    %   See also CACHE/IN, CACHE/OUT, CACHE/FIND, CACHE/COMPARE, CACHE/RESIZE,
    %       CACHE/CLEAR
    
    %   Andrew D. Horchler, adh9 @ case . edu, Created 11-29-13
    %   Revision: 1.0, 11-29-13
    
    
    properties (SetAccess=private)
        Size
        Used = 0;
    end
    
    properties (SetAccess=private,Hidden)
        numin
        numout
        inputs
        outputs
    end
    
    methods
        
        function obj = CACHE(Size,varargin)
            if nargin == 0
                Size = 1;
            elseif ~isscalar(Size) || ~isnumeric(Size) || ~isreal(Size) ...
                    || ~isfinite(Size) || Size < 1 || Size ~= floor(Size)
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
                IN(obj,varargin{:});
            end
        end
        
        function idx = IN(obj,varargin)
            idx = FIND(obj,varargin{:});
            if isempty(idx)
                if obj.Used < obj.Size
                    obj.Used = obj.Used+1;
                    obj.numin(obj.Used) = nargin-1;
                    obj.inputs{obj.Used} = varargin;
                else
                    obj.numin = [obj.numin(2:end);nargin-1];
                    [obj.inputs{:}] = deal(obj.inputs{2:end},varargin);
                end
            end
        end
        
        function [n,varargout] = OUT(obj,idx,varargin)
            if obj.Used == 0
                error('CACHE:OUT:NotInitialized',...
                   	  'Cache has not been initialized.');
            elseif nargin == 1 || isempty(idx)
             	idx = obj.Used;
        	elseif ~isscalar(idx) || ~isnumeric(idx) || ~isreal(idx) ...
                    || ~isfinite(idx) || idx < 1 || idx ~= floor(idx) || ...
                    idx > obj.Used
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
                    error('CACHE:OUT:TooManyOutputs',...
                          'Too many output arguments.');
                end
                v = 1:min(nargout-1,n);
                [varargout{v}] = obj.outputs{idx}{v};
            end
        end
        
        function idx = FIND(obj,varargin)
            idx = [];
            if ~isempty(varargin)
                for i = 1:obj.Used
                    if nargin-1 == obj.numin(i) ...
                            && isequaln(obj.inputs{i},varargin)
                        idx = i;
                        break;
                    end
                end
            end
        end
        
        function tf = COMPARE(obj,varargin)
            tf = ~isempty(FIND(obj,varargin{:}));
        end
        
        function obj = RESIZE(obj,Size)
            if nargin == 1
                Size = 1;
            elseif ~isscalar(Size) || ~isnumeric(Size) || ~isreal(Size) ...
                    || ~isfinite(Size) || Size < 0 || Size ~= floor(Size)
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
        end
        
    end
    
end