function net=shc_deletenetwork(net,i)
%SHC_DELETENETWORK  Remove a sub-network and its children.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 1-21-12
%   Revision: 1.0, 3-24-12


% Check network
shc_validatenetwork(net);

% Check network index
if ~validateindex(i) || ~isnumeric(i)
    error(  'SHCTools:shc_addnetwork:InvalidParentIndex',...
            'The index argument must be a real positive integer.');
end

if i == 1
    % Delete everything and return empty network
    net = struct();
    net.s = {};
else
    % Find index of first non-deleted subnetwork or end of network
    p = net.s{i}.parent;
    nnets = length(net.s);
    for j = i+1:nnets+1
        if j > nnets || net.s{j}.parent <= p
            break;
        end
    end

    % Remove deleted subnetwork(s)
    net.s = net.s([1:i-1 j:end]);

    % Reset network and re-initialize
    net = shc_initialize(net,'reset');
end