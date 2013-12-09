%% 3D Acoustic wave-equation
% 
% $$\left(\begin{array}{cc}\kappa^{-1}&0\\0&\rho \end{array}\right)\dot{\mathbf{w}} = \left(\begin{array}{cc}0&\nabla\cdot\\\nabla&0 \end{array}\right)\mathbf{w}$$
%


%% set parameters
%

% method, fourier or cheb
method = 'fourier';
% dimension
nd = 3;
% size of domain (m)
L = 1e3*ones(1,nd);
% # of gridpoints
N = 50*ones(1,nd);
% # of nodes for spectral method
Ns = 10*ones(1,nd);
% time interval
T = 1;
% medium parameters
rho0 = 1e3;
c0   = 1e3;

%% define matrices etc.
%

% Gradient operator
Grad = opGrad(Ns,L,method,true);

% spectral grid
xc = Grad.x{1};
yc = Grad.x{2};
zc = Grad.x{3};
[xxc,yyc,zzc] = ndgrid(xc,yc,zc);

% stiffness matrix, note that Div = -Grad' !
S    = [opZeros(prod(Ns)) Grad; -Grad' opZeros(nd*prod(Ns))];

% mass matrix
rho    = rho0*ones(Ns); 
kappa  = (c0^2*rho0)*ones(Ns);
M      = opDiag([kappa(:).^(-1);rho(:);rho(:);rho(:)]);

% regular grid
x = linspace(0,L(1),N(1));
y = linspace(0,L(2),N(2));
z = linspace(0,L(3),N(3));

% interpolation
A = opKron(opInt(zc,z,method),opInt(yc,y,method),opInt(xc,x,method));

% Initial conditions w = [p; ux; uy; uz]
w0        = zeros([Ns 4]);
w0(:,:,:,1) = exp(-1e3*(((xxc - mean(xc))/L(1)).^2 + ((yyc - mean(yc))/L(2)).^2 + ((zzc - mean(zc))/L(3)).^2));
w0        = w0(:)/max(abs(w0(:)));


%% solve ODE
tic
options  = odeset('Stats','on','OutputFcn',@odewbar);
[t,wsol] = ode23(@(t,w)(M\(S*w)),[0 T],w0,options) ;
toc

%% plot
wsol = reshape(wsol,[length(t) prod(Ns) 4]);
plotslice = @(u) slice(u,floor(N(1)/2) + 1,floor(N(2)/2) + 1,floor(N(3)/2) + 1);
for j=1:1:length(t);
    
    pj  = reshape(A*squeeze(wsol(j,:,1).'),N);
    u1j = reshape(A*squeeze(wsol(j,:,2).'),N);
    u2j = reshape(A*squeeze(wsol(j,:,3).'),N);
    u3j = reshape(A*squeeze(wsol(j,:,4).'),N);
    
    subplot(2,2,1);
    plotslice(pj);title('p');colormap(seiscol);
    subplot(2,2,2);
    plotslice(u1j);title('u_1');colormap(seiscol);
    subplot(2,2,3);
    plotslice(u2j);title('u_2');colormap(seiscol);
    subplot(2,2,4);
    plotslice(u3j);title('u_3');colormap(seiscol);
    drawnow;
    pause(.001);
end