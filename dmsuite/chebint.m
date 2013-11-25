function p = chebint(fk, x, mode)

%  The function p = chebint(fk, x) computes the polynomial interpolant
%  of the data (xk, fk), where xk are the Chebyshev nodes.  
%  Two or more data points are assumed.
%
%  Input:
%  fk:  Vector of y-coordinates of data, at Chebyshev points 
%       x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  x:   Vector of x-values where polynomial interpolant is to be evaluated.
%  mode:     1:forward, -1:adjoint
%
%  Output:
%  p:    Vector of interpolated values.
%
%  The code implements the barycentric formula; see page 252 in
%  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)

%  J.A.C. Weideman, S.C. Reddy 1998
%
% Modification by Tristan van Leeuwen, November 2013:
% - added flag for forward and adjoint mode
% - works on input matrix, doing the same for each column.

x = x(:);                    % Make sure data are column vectors.

[N, Nc]  = size(fk);
M        = length(x);

% Compute Chebyshev points.
xk = sin(pi*[N-1:-2:1-N]'/(2*(N-1)));    

% w = weights for Chebyshev formula
w    = ones(N,1).*(-1).^[0:N-1]';         
w(1) = w(1)/2; w(N) = w(N)/2;

% Compute quantities x-x(k) and their reciprocals.
D = x(:,ones(1,N)) - xk(:,ones(1,M))'; 
D(D==0) = eps;
D = 1./D;

% Evaluate interpolant as matrix-vector products.
a = D*w;
if mode == 1
    fk = bsxfun(@times,fk,w);
    p  = D*fk;    
    p  = bsxfun(@times,p,1./a);
else
    fk = bsxfun(@times,fk,1./a);
    p  = D'*fk;
    p  = bsxfun(@times,p,w);
end
