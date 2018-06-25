%% 3D Maxwell's equation
% 
% $$\left(\begin{array}{cc}\epsilon&0\\0&\mu \end{array}\right)\dot{\mathbf{w}} + \left(\begin{array}{cc}0&\nabla\times\\-\nabla\times&0 \end{array}\right)\mathbf{w} = 0$$
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
N = 30*ones(1,nd);
% # of nodes for spectral method
Ns = 10*ones(1,nd);
% time interval
T = 1e-6;
% medium parameters
eps0 = 8.85e-12; % di-electric permitivity 
mu0  = 1.25e-6;  % magnetic permeability

%% define matrices etc.
%
store = true;
% Grad
Grad = opGrad(Ns,L,method,store);
% Curl
Curl = opCurl(Ns,L,method,store);
%

% spectral grid
xc = Grad.x{1};
yc = Grad.x{2};
zc = Grad.x{3};
[xxc,yyc,zzc] = ndgrid(xc,yc,zc);

% stiffness matrix
S    = [opZeros(nd*prod(Ns)) Curl; -Curl opZeros(nd*prod(Ns))];

% mass matrix
eps    = eps0*ones(Ns); 
mu     = (mu0)*ones(Ns);
M      = opDiag([eps(:);eps(:);eps(:);mu(:);mu(:);mu(:)]);

% regular grid
x = linspace(0,L(1),N(1));
y = linspace(0,L(2),N(2));
z = linspace(0,L(3),N(3));

% interpolation
A = opKron(opInt(zc,z,method),opInt(yc,y,method),opInt(xc,x,method));

% Initial conditions w = [ex; ey; ex; hx; hy; hz]
w0        = zeros([Ns 6]);
w0(:,:,:,1) = (yyc - mean(yc)).*exp(-1e3*(((xxc - mean(xc))/L(1)).^2 + ((yyc - mean(yc))/L(2)).^2 + ((zzc - mean(zc))/L(3)).^2));
w0(:,:,:,2) = -(xxc - mean(xc)).*exp(-1e3*(((xxc - mean(xc))/L(1)).^2 + ((yyc - mean(yc))/L(2)).^2 + ((zzc - mean(zc))/L(3)).^2));
w0          = w0(:)/max(abs(w0(:)));


%% solve ODE
tic
options  = odeset('Stats','on','OutputFcn',@odewbar);
[t,wsol] = ode23(@(t,w)(M\(S*w)),[0 T],w0,options) ;
toc

%% plot
wsol      = reshape(wsol,[length(t) prod(Ns) 6]);
plotslice = @(u) slice(u,floor(N(1)/2) + 1,floor(N(2)/2) + 1,floor(N(3)/2) + 1);
for j=1:10:length(t);
    
    exj = reshape(A*squeeze(wsol(j,:,1).'),N);
    eyj = reshape(A*squeeze(wsol(j,:,2).'),N);
    ezj = reshape(A*squeeze(wsol(j,:,3).'),N);
    
    hxj = reshape(A*squeeze(wsol(j,:,1).'),N);
    hyj = reshape(A*squeeze(wsol(j,:,2).'),N);
    hzj = reshape(A*squeeze(wsol(j,:,3).'),N);
    
    subplot(2,3,1);
    plotslice(exj);title('E_x');colormap(seiscol);
    subplot(2,3,2);
    plotslice(eyj);title('E_y');colormap(seiscol);
    subplot(2,3,3);
    plotslice(ezj);title('E_z');colormap(seiscol);

    subplot(2,3,4);
    plotslice(hxj);title('H_x');colormap(seiscol);
    subplot(2,3,5);
    plotslice(hyj);title('H_y');colormap(seiscol);
    subplot(2,3,6);
    plotslice(hzj);title('H_z');colormap(seiscol);
    
    drawnow;
    pause(1);
end