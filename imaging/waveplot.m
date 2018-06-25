function waveplot(wsol,params)

%% set parameters
method = params.method;
nd     = params.nd;
L      = params.L;
N      = params.N;
Ns     = params.Ns;
T      = params.T;
dt     = params.dt;
t      = [0:dt:T];
%% define operators
Grad = opGrad(Ns,L,method,false);

% spectral grid
xc = Grad.x{1};
yc = Grad.x{2};

% regular grid
x = linspace(0,L(1),N(1));
y = linspace(0,L(2),N(2));

% interpolation
A = opKron(opInt(yc,y,method),opInt(xc,x,method));

%% plot
wsol = reshape(wsol,[prod(Ns) 3 length(t)]);
for j=1:ceil(length(t)/100):length(t);
    pj = reshape(A*squeeze(wsol(:,1,j)),N);
    imagesc(y,x,pj,[-1 1]);title('p');colormap(seiscol);axis equal tight;
    drawnow;
    pause(.001);
end