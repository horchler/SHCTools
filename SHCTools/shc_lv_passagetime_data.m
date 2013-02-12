function varargout=shc_lv_passagetime_data(net,delta,t,a)
%SHC_LV_PASSAGETIME_DATA  
%
%   TAU = SHC_LV_PASSAGETIME_DATA(NET,DELTA,T,A)
%   [TP,TT] = SHC_LV_PASSAGETIME_DATA(NET,DELTA,T,A)
%   [TAU,TP,TT] = SHC_LV_PASSAGETIME_DATA(NET,DELTA,T,A)
%   [...,TI] = SHC_LV_PASSAGETIME_DATA(NET,DELTA,T,A)

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-2-12
%   Revision: 1.0, 2-12-13


if nargout > 4
    error('SHCTools:shc_lv_passagetime_data:TooManyOutputs',...
          'Too many output arguments.');
end

% Check network
if ~isstruct(net) || ~isfield(net,'rho')
    error('SHCTools:shc_lv_passagetime_data:NetworkStructOrRhoInvalid',...
          'Input must be a valid SHC network structure.');
end
n = net.size;

% Check Delta
if ~isvector(delta) || isempty(delta) || ~isfloat(delta)
    error('SHCTools:shc_lv_passagetime_data:DeltaInvalid',...
         ['The neighborhood size, Delta, must be a non-empty floating-point '...
          'vector.']);
end
if ~isscalar(delta) && length(delta) ~= n
    error('SHCTools:shc_lv_passagetime_data:DeltaDimensionMismatch',...
         ['If the neighborhood size, Delta, is a vector, it must have the '...
          'same length as the network dimension.']);
end
if ~isreal(delta) || ~all(isfinite(delta)) || any(delta <= 0)
    error('SHCTools:shc_lv_passagetime_data:DeltaNonFiniteReal',...
         ['The neighborhood size, Delta, must be a positive finite real '...
          'floating-point vector.']);
end
delta = delta(:);

% Handle inputs for A and T
if ~shc_ismatrix(a) || ~isreal(a) || ~isfloat(a) || ~all(isfinite(a(:)))
    error('SHCTools:shc_lv_passagetime_data:InvalidA',...
         ['The data, A, must be a finite real matrix of floating point '...
          'values.']);
end
[M,N] = size(a);
if N ~= n || M < 2
    error('SHCTools:shc_lv_passagetime_data:DimensionMismatchA',...
         ['The data matrix, A, must have at least two rows and number of '...
          'columns must equal the network size.']);
end

if ~isvector(t) || ~isreal(t) || ~isfloat(t) || ~all(isfinite(t))
    error('SHCTools:shc_lv_passagetime_data:InvalidT',...
         ['The time data, T, must be a finite real vector of floating point '...
          'values.']);
end
if length(t) ~= M
    error('SHCTools:shc_lv_passagetime_data:DimensionMismatchT',...
         ['The time vector, T, must have the same length as the number of '...
          'rows of the data matrix, A.']);
end

dagtd = diff(a > delta)';
[j1,i1] = find([false(n,1) (dagtd < 0)]);
[j2,i2] = find([false(n,1) (dagtd > 0)]);
if i1(1) > i2(1)
    i2 = i2(2:end);
    j2 = j2(2:end);
end
lt = min(length(i1)-1,length(i2));
if lt > 0
    if nargout ~= 2
        tau = num2cell(NaN(lt,n),1);
    end
    if nargout ~= 1
        tp = num2cell(NaN(lt,n),1);
        tt = num2cell(NaN(lt,n),1);
    end
    if nargout == 4
        ti = num2cell(NaN(lt,n),1);
    end
    ic(1,n) = 0;
    ti1 = t(i1(1)-1)+(t(i1(1))-t(i1(1)-1))*(delta-a(i1(1)...
        -1,j1(1)))/(a(i1(1),j1(1))-a(i1(1)-1,j1(1)));
    for i = 1:lt
        j1i = j1(i);
        ic(j1i) = ic(j1i)+1;
        icj = ic(j1i);
        
        i1i2 = i1(i+1);
        j1i2 = j1(i+1);
        ti12 = t(i1i2-1)+(t(i1i2)-t(i1i2-1))*(delta-a(i1i2...
            -1,j1i2))/(a(i1i2,j1i2)-a(i1i2-1,j1i2));
        if nargout ~= 1
            i2i = i2(i);
            j2i = j2(i);
            ti2 = t(i2i-1)+(t(i2i)-t(i2i-1))*(delta-a(i2i-1,j2i))/(a(i2i,j2i)...
                -a(i2i-1,j2i));
            tp{j1i}(icj) = ti2-ti1;
            tt{j1i}(icj) = ti12-ti2;
            if nargout == 4
                ti{j1i}(icj) = ti1;
            end
        end
        if nargout ~= 2
            tau{j1i}(icj) = ti12-ti1;
        end
        ti1 = ti12;
    end
    
    % Handle variable output and remove any remaining NaNs
    if nargout == 1
        varargout{1} = cellfun(@(x)x(~isnan(x),~all(isnan(x))),tau,...
            'UniformOutput',false);
    elseif nargout == 2
        varargout{1} = cellfun(@(x)x(~isnan(x),~all(isnan(x))),tp,...
            'UniformOutput',false);
        varargout{2} = cellfun(@(x)x(~isnan(x),~all(isnan(x))),tt,...
            'UniformOutput',false);
    elseif nargout > 2
        varargout{1} = cellfun(@(x)x(~isnan(x),~all(isnan(x))),tau,...
            'UniformOutput',false);
        varargout{2} = cellfun(@(x)x(~isnan(x),~all(isnan(x))),tp,...
            'UniformOutput',false);
        varargout{3} = cellfun(@(x)x(~isnan(x),~all(isnan(x))),tt,...
            'UniformOutput',false);
        if nargout == 4
            varargout{4} = cellfun(@(x)x(~isnan(x),~all(isnan(x))),ti,...
                'UniformOutput',false);
        end
    end
else
    if nargout == 1
        varargout{1} = cell(1,n);
    elseif nargout == 2
        varargout{1} = cell(1,n);
        varargout{2} = cell(1,n);
    elseif nargout > 2
        varargout{1} = cell(1,n);
        varargout{2} = cell(1,n);
        varargout{3} = cell(1,n);
        if nargout == 4
            varargout{4} = cell(1,n);
        end
    end
end