WaveMat
=======

Spectral methods for wave propagation in Matlab. The code is designed to 
solve coupled ODE's of the form

\ M\mathbf{w}(t) + S\mathbf{w}(t) = \mathbf{f}(t) \

where the spatial differenatiation in \S\ is done via Spectral methods. The code is basically a wrapper
for the "Differentiation Matrix Suite" (http://dip.sun.ac.za/~weideman/research/differ.html) using 
the SPOT framework (http://www.cs.ubc.ca/labs/scl/spot/). This allows one to effortlessly define
Mass and Stiffness matrices, which can then be passed to any ODE solved such as Matlab's `odeXY`.

Some examples are included in `/examples`.
