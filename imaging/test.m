%% 2D Acoustic wave-equation
% 
% $$\left(\begin{array}{cc}\kappa^{-1}&0\\0&\rho \end{array}\right)\dot{\mathbf{w}} = \left(\begin{array}{cc}0&\nabla\cdot\\\nabla&0 \end{array}\right)\mathbf{w}$$
%
% 

%% set parameters
%

% method, fourier or cheb
params.method = 'cheb';
% dimension
params.nd = 2;
% size of domain (m)
params.L = 1e3*ones(1,nd);
% # of gridpoints
params.N = 100*ones(1,nd);
% # of nodes for spectral method
params.Ns = 50*ones(1,nd);
% time interval
params.T = .5;
params.dt = 1e-3;
% damping
params.beta = 1e3;
params.Npml = 25*ones(1,nd);
% medium parameters
c = 1e3*ones(params.Ns); % velocity in m/s

%% initial condition
w0 = PointSource([500 50],params);

%% solve ODE
wsol =  MFull(c,w0,1,params);

%% plot
waveplot(wsol,params)