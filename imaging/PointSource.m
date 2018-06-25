function w0 = PointSource(p,params)

%% set parameters
method = params.method;
nd     = params.nd;
L      = params.L;
N      = params.N;
Ns     = params.Ns;
T      = params.T;
dt     = params.dt;

%% define operators
Grad = opGrad(Ns,L,method,false);

% spectral grid
xc = Grad.x{1};
yc = Grad.x{2};
[xxc,yyc] = meshgrid(xc,yc);
% regular grid
%x = linspace(0,L(1),N(1));
%y = linspace(0,L(2),N(2));

% interpolation
%A = opKron(opInt(yc,y,method),opInt(xc,x,method));

%%
w0 = zeros([Ns 3]);
w0(:,:,1) = exp(-1e3*(((xxc - p(1))/L(1)).^2 + ((yyc - p(2))/L(2)).^2));
w0 = w0(:)/max(abs(w0(:)));
