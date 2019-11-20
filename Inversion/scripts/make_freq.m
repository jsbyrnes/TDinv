function f=make_freq(n,dt)
%   make_freq     make frequency vector for Fourier Transforms
% USAGE: f=make_freq(n,dt);
%  construct the frequency vector that is consistent with the result of an FFT.
%  dt is the time-domain sample interval and n is the number of points.
%
%  See also FT and IFT.

% K. Creager  kcc@geophys.washington.edu   12/30/93

if rem(n,2)==1
  f=[0:(n-1)/2  -(n-1)/2:-1]'/(n*dt);   % frequency vector for odd n
else
  f=[0:n/2        -n/2+1:-1]'/(n*dt);   % frequency vector for even n
end
