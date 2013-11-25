function Dmf = chebdifft(f,M,mode)

% The function Dmf = chebdifft(f,M) computes the M'th
% approximate Chebyshev derivatives of the data vector y.
% A Fast Fourier Transform is used compute the Chebyshev cofficients
% of the data vector. A recursion formula is used to compute the
% Chebyshev coefficients for each derivative. A FFT is then used again
% to compute the derivatives in physical space.
% 
%  Input:
%  f:        Vector containing function values at the Chebyshev points
%            x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  mode:     1:forward, -1:adjoint
%  M:        Derivative required (positive integer)
%
%  Output:
%  ym:       Vector containing approximate M'th derivative
%            
%  J.A.C. Weideman, S.C. Reddy 2000.  
%
% Modification by Tristan van Leeuwen, November 2013:
% - added flag for forward and adjoint mode
% - works on input matrix, doing the same for each column.

N  = size(f,1);
Nc = size(f,2);

if mode == 1
    f = [f; flipud(f(2:N-1,:))];
    a = fft(f);
    a = (a(1:N,:) + [zeros(1,Nc);flipud(a(N+1:end,:));zeros(1,Nc)])/(2*(N-1));
    
    A = spdiags([[ones(N-1,1);1] -[.5;.5;.5;ones(N-3,1)]],[0 2],N,N);
    B = spdiags([1;1;2*[2:N]'],1,N,N);
    for ell=1:M
        a = A\(B*a);
    end
    
    a    = [2*a(1,:); a(2:N-1,:); 2*a(N,:)];
    a    = [a; flipud(a(2:N-1,:))];
    Dmf  = .5*fft(a);
    Dmf  = (Dmf(1:N,:) + [zeros(1,Nc);flipud(Dmf(N+1:end,:));zeros(1,Nc)]);
    Dmf  = [Dmf(1,:); .5*Dmf(2:N-1,:); Dmf(N,:)];
else
   f = [f(1,:); .5*f(2:N-1,:); f(N,:)];
   f = [f; flipud(f(2:N-1,:))];
   a = .5*fft(f);
   a = a(1:N,:) + [zeros(1,Nc);flipud(a(N+1:end,:));zeros(1,Nc)];
   a = [2*a(1,:); a(2:N-1,:); 2*a(N,:)];
    
    A = spdiags([[ones(N-1,1);1] -[.5;.5;.5;ones(N-3,1)]],[0 2],N,N);
    B = spdiags([1;1;2*[2:N]'],1,N,N);
    for ell=1:M
        a = B'*(A'\a);
    end
    
    a    = [a; flipud(a(2:N-1,:))]/(2*(N-1));
    Dmf  = fft(a);
    Dmf  = Dmf(1:N,:) + [zeros(1,Nc);flipud(Dmf(N+1:end,:));zeros(1,Nc)];

end

Dmf = real(Dmf);
