%% 2D Acoustic wave-equation with PML
% 
% $$\left(\begin{array}{cccc}\kappa^{-1}&0&0&0\\0&\kappa^{-1}&0&0\\0&0&\rho &0\\0&0&0&\rho \end{array}\right)\dot{\mathbf{w}} = \left(\begin{array}{cccc}0&0&\partial_x&0\\0&0&0&\partial_y \\ \partial_x&\partial_x&0&0\\ \partial_y&\partial_y&0&0 \end{array}\right)\mathbf{w} - \left(\begin{array}{cccc}\sigma_x&0&0&0\\ 0&\sigma_y&0&0\\0&0&\sigma_x &0\\0&0&0&\sigma_y \end{array}\right)\mathbf{w}$$
%
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
% # of PML points
Npml = 10*ones(1,nd);
% PML strength
beta = 100;
% time interval
T = 1;
% medium parameters
rho0 = 1e3; % density kg/m^3
c0   = 1e3; % velocity in m/s

%% define matrices etc.
%

% Gradient operators
Dx = opDel(Ns,L,1,method,true);
Dy = opDel(Ns,L,2,method,true);

%
O = opZeros(prod(Ns));

% spectral grid
xc = Dx.x{1};
yc = Dy.x{2};
[xxc,yyc] = ndgrid(xc,yc);

% stiffness matrix, note that Div = -Grad' !
S    = [O O Dx O;O O O Dy;-Dx' -Dx' O O;-Dy' -Dy' O O];

% mass matrix
rho    = rho0*ones(Ns); 
kappa  = (c0^2*rho0)*ones(Ns);
M      = opDiag([kappa(:).^(-1);kappa(:).^(-1);rho(:);rho(:)]);

% PML part
sigmax = [beta*linspace(1,0,Npml(1)).^2 zeros(1,Ns(1)-2*Npml(1)) beta*linspace(0,1,Npml(1)).^2]'*ones(1,Ns(2));
sigmay = ones(Ns(1),1)*[beta*linspace(1,0,Npml(2)).^2 zeros(1,Ns(2)-2*Npml(2)) beta*linspace(0,1,Npml(2)).^2];

B = opDiag([sigmax(:);sigmay(:);sigmax(:);sigmay(:)]);


% regular grid
x = linspace(0,L(1),N(1));
y = linspace(0,L(2),N(2));

% interpolation
A = opKron(opInt(yc,y,method),opInt(xc,x,method));

% Initial conditions, w = [p; ux ; uy]
w0        = zeros([Ns 4]);
w0(:,:,1) = exp(-1e3*(((xxc - mean(xc))/L(1)).^2 + ((yyc - mean(yc))/L(2)).^2));
w0(:,:,2) = exp(-1e3*(((xxc - mean(xc))/L(1)).^2 + ((yyc - mean(yc))/L(2)).^2));
w0        = w0(:)/max(abs(w0(:)));


%% solve ODE
tic
options  = odeset('Stats','on','OutputFcn',@odewbar);
[t,wsol] = ode23(@(t,w)(M\(S*w - M*B*w)),[0 T],w0,options) ;
toc

%% plot
wsol = reshape(wsol,[length(t) prod(Ns) 4]);

for j=1:1:length(t)
    
    pj = reshape(A*squeeze(wsol(j,:,1).' + wsol(j,:,2).'),N);
    uj = reshape(A*squeeze(wsol(j,:,3).'),N);
    vj = reshape(A*squeeze(wsol(j,:,4).'),N);
    
    subplot(1,3,1);
    imagesc(y,x,pj,[-1 1]*1e-1);title('p');colormap(seiscol);axis equal tight;
    subplot(1,3,2);
    imagesc(y,x,uj,[-1 1]*1e-7);title('u');colormap(seiscol);axis equal tight;
    subplot(1,3,3);
    imagesc(y,x,vj,[-1 1]*1e-7);title('v');colormap(seiscol);axis equal tight;  
    drawnow;
    pause(.001);
end