function wsol = MFull(c,w0,dir,params)
% propagator based on full wave equation
%
% use:
%    u = MFull(c,w0,dir,params)
%
% input:
%   c  - soundspeed prod(N) x 1 array
%   w0 - initial condition prod(Ns) x 3 array
%   params.method - differenatiation {cheb, fourier, fd};
%   params.nd     - number of spatial dimensions;
%   params.L      - size of domain;
%   params.N      - number of gridpoints for c;
%   params.Ns     - number of gridpoints for differentation;
%   params.T      - time;
%   params.dt     - timestep;
%
% output:
%   u - u(t,x)
%
%% set parameters
method = params.method;
nd     = params.nd;
L      = params.L;
N      = params.N;
Ns     = params.Ns;
T      = params.T;
dt     = params.dt;
beta   = params.beta;
Npml   = params.Npml;

%% define operators
% Gradient operator
Grad = opGrad(Ns,L,method,false);

% stiffness matrix, note that Div = -Grad' !
S    = [opZeros(prod(Ns)) -Grad'; Grad opZeros(nd*prod(Ns))];

% mass matrix
rho0   = 1;
rho    = rho0*ones(Ns); 
kappa  = (c.^2.*rho0);
M      = opDiag([kappa(:).^(-1);rho(:);rho(:)]);

% spectral grid
xc = Grad.x{1};
yc = Grad.x{2};

% regular grid
x = linspace(0,L(1),N(1));
y = linspace(0,L(2),N(2));

% interpolation
A = opKron(opInt(yc,y,method),opInt(xc,x,method));

% Damping part
sigmax = [beta*linspace(1,0,Npml(1)).^2 zeros(1,Ns(1)-2*Npml(1)) beta*linspace(0,1,Npml(1)).^2]'*ones(1,Ns(2));
sigmay = ones(Ns(1),1)*[beta*linspace(1,0,Npml(2)).^2 zeros(1,Ns(2)-2*Npml(2)) beta*linspace(0,1,Npml(2)).^2];
sigma  = sigmax + sigmay;

B = opDiag([sigma(:);sigma(:);sigma(:)]);

%% timestepping
tic
%[t,wsol] = ode23(@(t,w)(M\(S*w)),[0:dt:T],w0) ;
t = 0:dt:T;
wsol = zeros([size(w0) length(t)]);
wsol(:,1) = w0;
for k = 2:length(t)
    wsol(:,k) = wsol(:,k-1) + dt*(M\(S*wsol(:,k-1)) - B*wsol(:,k-1));
end
toc

%% reshape
wsol = reshape(wsol,[length(t) prod(Ns) 3]);
