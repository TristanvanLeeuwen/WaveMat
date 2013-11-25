    function Dmf = fourdifft(f,m,mode)

% Dmf = fourdifft(f,m) computes the m-th derivative of the function
% f(x) using the Fourier differentiation process.   The function 
% is assumed to be 2pi-periodic and the input data values f should 
% correspond to samples of the function at N equispaced points on [0, 2pi).
% The Fast Fourier Transform is used.
% 
%  Input:
%  f:      Vector of samples of f(x) at x = 0, 2pi/N, 4pi/N, ... , (N-1)2pi/N
%  m:      Derivative required (non-negative integer)
%  mode:   1:forward, -1:adjoint
%  Output:
%  Dmf:     m-th derivative of f
%
%  S.C. Reddy, J.A.C. Weideman 2000. Corrected for MATLAB R13 
%  by JACW, April 2003.
%
% Modification by Tristan van Leeuwen, November 2013:
% - added flag for forward and adjoint mode
% - works on input matrix, doing the same for each column.


%f = f(:);                       % Make sure f is a column vector
N  = size(f,1);
N1 = floor((N-1)/2);           % Set up wavenumbers
N2 = (-N/2)*rem(m+1,2)*ones(rem(N+1,2));
wave = [(0:N1)  N2 (-N1:-1)]';

if mode==1
    Dmf = ifft(bsxfun(@times,(+1i*wave).^m,fft(f)));   % Transform to Fourier space, take deriv, and return to physical space.
else
    Dmf = ifft(bsxfun(@times,(-1i*wave).^m,fft(f)));  % adjoint mode
end

if max(abs(imag(f))) == 0; 
    Dmf = real(Dmf);  % Real data in, real derivative out
end 
