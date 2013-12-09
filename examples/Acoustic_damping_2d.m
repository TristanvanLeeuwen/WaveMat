%% 2D Acoustic wave-equation with Damping boundary
% 
% $$\left(\begin{array}{cc}\kappa^{-1}&0\\0&\rho \end{array}\right)\dot{\mathbf{w}} = \left(\begin{array}{cc}0&\nabla\cdot\\\nabla&0 \end{array}\right)\mathbf{w} - \left(\begin{array}{cc}\sigma&0\\0&\sigma \end{array}\right)\mathbf{w}$$
%

%% set parameters
%

% method, fourier or cheb
method = 'fourier';
% dimension
nd = 2;
% size of domain (m)
L = 1e3*ones(1,nd);
% # of gridpoints
N = 100*ones(1,nd);
% # of nodes for spectral method
Ns = 50*ones(1,nd);
% # of damping points
Npml = 10*ones(1,nd);
% Damping strength
beta = 100;
% time interval
T = 1;
% medium parameters
rho0 = 1e3; % density kg/m^3
c0   = 1e3; % velocity in m/s

%% define matrices etc.
%

% Gradient operator
Grad = opGrad(Ns,L,method,true);

% spectral grid
xc = Grad.x{1};
yc = Grad.x{2};
[xxc,yyc] = ndgrid(xc,yc);

% stiffness matrix, note that Div = -Grad' !
S    = [opZeros(prod(Ns)) Grad; -Grad' opZeros(nd*prod(Ns))];

% mass matrix
rho    = rho0*ones(Ns); 
kappa  = (c0^2*rho0)*ones(Ns);
M      = opDiag([kappa(:).^(-1);rho(:);rho(:)]);

% Damping part
sigmax = [beta*linspace(1,0,Npml(1)).^2 zeros(1,Ns(1)-2*Npml(1)) beta*linspace(0,1,Npml(1)).^2]'*ones(1,Ns(2));
sigmay = ones(Ns(1),1)*[beta*linspace(1,0,Npml(2)).^2 zeros(1,Ns(2)-2*Npml(2)) beta*linspace(0,1,Npml(2)).^2];
sigma  = sigmax + sigmay;

B = opDiag([sigma(:);sigma(:);sigma(:)]);


% regular grid
x = linspace(0,L(1),N(1));
y = linspace(0,L(2),N(2));

% interpolation
A = opKron(opInt(yc,y,method),opInt(xc,x,method));

% Initial conditions, w = [p; ux ; uy]
w0        = zeros([Ns 3]);
w0(:,:,1) = exp(-1e3*(((xxc - mean(xc))/L(1)).^2 + ((yyc - mean(yc))/L(2)).^2));
w0        = w0(:)/max(abs(w0(:)));


%% solve ODE
tic
options  = odeset('Stats','on','OutputFcn',@odewbar);
[t,wsol] = ode23(@(t,w)(M\(S*w - M*B*w)),[0 T],w0,options) ;
toc

%% plot
wsol = reshape(wsol,[length(t) prod(Ns) 3]);

for j=1:1:length(t);
    
    pj = reshape(A*squeeze(wsol(j,:,1).'),N);
    uj = reshape(A*squeeze(wsol(j,:,2).'),N);
    vj = reshape(A*squeeze(wsol(j,:,3).'),N);
    
    subplot(1,3,1);
    imagesc(y,x,pj,[-1 1]*1e-1);title('p');colormap(seiscol);axis equal tight;
    subplot(1,3,2);
    imagesc(y,x,uj,[-1 1]*1e-7);title('u');colormap(seiscol);axis equal tight;
    subplot(1,3,3);
    imagesc(y,x,vj,[-1 1]*1e-7);title('v');colormap(seiscol);axis equal tight;  
    drawnow;
    pause(.001);
end