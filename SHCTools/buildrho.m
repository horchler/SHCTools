function net=buildrho(net)
%BUILDRHO  Create RHO matrix from SHC network structure.
%
%

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-21-10
%   Revision: 1.0, 6-30-12


% Convert from XML file to structure if necessary
if ischar(net)
    net = loadnet(net);
    net = shc_initialize(net);
end

% Validate network structure
shc_validatenetwork(net);

m = net.size;
n = net.s{1}.size;
z(n,1) = 0;
zz = zeros(m-n,1);

% Initialize parameter vectors
net.alpha = [z+net.s{1}.alpha;zz];
net.beta = [z+net.s{1}.beta;zz];
net.gamma = [z+net.s{1}.gamma;zz];
net.delta = [z+net.s{1}.delta;zz];

% Initialize T matrix
net.T = [net.s{1}.T false(n,m-n);false(m-n,m)];

% Build T matrix
for i = 2:length(net.s)
    m = n;
    n = net.s{i}.size;
    x = net.s{i}.index;
    if isempty(zz)
        z = 0;
    else
        z = zz(1:n);
    end
    
    n = n+m-1;
    net.alpha(x:n) = [z+net.s{i}.alpha;net.alpha(x+1:m)];
    net.beta(x:n) = [z+net.s{i}.beta;net.beta(x+1:m)];
    net.gamma(x:n) = [z+net.s{i}.gamma;net.gamma(x+1:m)];
    net.delta(x:n) = [z+net.s{i}.delta;net.delta(x+1:m)];
    
    if strncmp(net.s{i}.type,'cluster',2)
        q = z+x;
        
        qx1 = net.T(q,1:x-1);
        qx2 = net.T(q,x+1:m);
        qy1 = net.T(1:x-1,q);
        qy2 = net.T(x+1:m,q);
    else
        q1 = false(n-m,x-1);
        q2 = false(n-m,m-x);
        
        if net.s{i}.direction == 1
            qx1 = [net.T(x,1:x-1);q1];
            qx2 = [net.T(x,x+1:m);q2];
            qy1 = [q1' net.T(1:x-1,x)];
            qy2 = [q2' net.T(x+1:m,x)];
        else
            qx1 = [q1;net.T(x,1:x-1)];
            qx2 = [q2;net.T(x,x+1:m)];
            qy1 = [net.T(1:x-1,x) q1'];
            qy2 = [net.T(x+1:m,x) q2'];
        end
    end
    
    % Build T matrix
    net.T(1:x-1,x:n) = [qy1 net.T(1:x-1,x+1:m)];
    net.T(x:n,1:n) = [qx1                   net.s{i}.T	qx2;
                      net.T(x+1:m,1:x-1)	qy2         net.T(x+1:m,x+1:m)];
end

% convert T to rho
q = ones(1,n);
net.rho = ~net.T.*net.gamma(:,q)+net.T.*net.delta(:,q);
net.rho(1:n+1:end) = net.alpha./net.beta;
if nargout == 0
    net = net.rho;
end