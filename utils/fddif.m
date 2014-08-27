function [x,D] = fddif(N,l)
% generate finite-difference differentiation matrix on interval [0,1]

h = 1/(N - 1);
x = linspace(0,1,N);
D = (1/h)*spdiags(ones(N,1)*[-1 1],[0 1],N,N);