function [t,wout] = fEuler(fh,TSPAN,w,dt,dsave)
% forward Euler for solving y' = f(t,y)
%
% use:
%   [t,y] = fEuler(f,TSPAN,y0,dt,dsave)
%
% input:
%   f     - function handle f(t,y) 
%   TSPAN - [Tmin Tmax]
%   y0    - initial condition
%   dt    - time-step
%   dsave -  
%
% output
%   t - time vector
%   y - solution at each time-step
%

t  = TSPAN(1):dt:TSPAN(2);
nt = length(t);
nsave = ceil(nt/dsave);
wout  = zeros(length(w),nsave);
wout(:,1) = w;
l = 2;
h = waitbar(0,'Time-stepping...');
for k = 1:nt-1
    w = w + dt*fh(t(k),w);
    if ~mod(k,dsave)
        wout(:,l) = w;
        l = l + 1;
    end
    waitbar(k/nt,h);
end
close(h);

t    = t(1:dsave:end);
wout = transpose(wout);