function net=buildrho(net)
%BUILDRHO  Create RHO matrix from SHC network structure.
%
%   See also:
%       SHC_CREATENETWORK, SHC_INITIALIZE, SHC_CREATE, SHC__LV_PARAMS,
%       SHC_LV_SYMPARAMS

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-21-10
%   Revision: 1.2, 5-4-13


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

% Initialize parameter vectors, gamma matrix
net.alpha = [z+net.s{1}.alpha;zz];
net.beta = [z+net.s{1}.beta;zz];
if isa(net.s{1}.gamma,'sym')
    net.gamma = sym(zeros(m,m));
else
    net.gamma(m,m) = 0;
end
if isscalar(net.s{1}.gamma)
    gam = double(~net.s{1}.T);
    gam(1:n+1:end) = 0;
    net.gamma(1:n,1:n) = net.s{1}.gamma.*gam;
elseif isvector(net.s{1}.gamma)
  	gam = double(~net.s{1}.T);
    gam(1:n+1:end) = 0;
    net.gamma(1:n,1:n) = net.s{1}.gamma(:,ones(1,n)).*gam;
else
    net.gamma(1:n,1:n) = net.s{1}.gamma;
end
net.delta = [z+net.s{1}.delta;zz];
if isfield(net.s{1},'nu')
    nu = [z+net.s{1}.nu;zz];
    isNu = true;
else
    isNu = false;
end

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
    
    pgam = net.gamma(x+1:m,x+1:m);
    if isa(net.s{i}.gamma,'sym') && ~isa(net.gamma,'sym')
        net.gamma = sym(net.gamma);
    end
    net.gamma(x:n,x:n) = 0;
    if isscalar(net.s{i}.gamma)
        gam = net.s{i}.gamma.*double(~net.s{i}.T);
    elseif isvector(net.s{i}.gamma)
        gam = double(~net.s{i}.T);
        gam(1:n-m+2:end) = 0;
        gam = net.s{i}.gamma(:,ones(1,n-m+1)).*gam;
    else
        gam = net.s{i}.gamma;
    end
    net.gamma(x+n+1:n,x+n+1:n) = pgam;
    net.gamma(x:x+n-m,x:x+n-m) = gam;
    
    net.delta(x:n) = [z+net.s{i}.delta;net.delta(x+1:m)];
    if isNu && isfield(net.s{i},'nu')
        nu(x:n) = [z+net.s{i}.nu;nu(x+1:m)];
    else
        isNu = false;
    end
    
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
if isNu
    net.nu = nu';
end

% convert T to rho
gam2 = ~net.T.*net.gamma;
net.rho = gam2+net.T.*net.delta(:,ones(1,m));
net.rho(1:m+1:end) = net.alpha./net.beta;

% convert Gamma matrix to column vector if possible
gamv = true;
for i = m:-1:1
    gami = gam2(i,gam2(i,:) ~= 0);
    gam3(i) = gami(1);
    if any(gami(1) ~= gami)
        gamv = false;
        break;
    end
end
if gamv
    net.gamma = gam3(:);
end

if nargout == 0
    net = net.rho;
end