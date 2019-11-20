% Input file for mat_disperse.m
clear, close all

% Define vectors or layer thickness, mass density, shear wave velocity
% and compression wave velocity
thk = ones(40,1);
dns = 1.8*ones(41,1);
vs = 4000*ones(41,1);
vp = 6000*ones(41,1);

vs(1:10) = 20;
vp(1:10) = 40;

% Define a vector of frequencies (in Hz)
freq = linspace(0.05,10,20);

% Define a vector of offsets from the source
offsets = 1e9;%linspace(5,100,20);

% Call mat_disperse.m to solve the eigenvalue problem and calculate phase
% velocities, displacement-stress functions, and surface wave displacements
tic
[vr,z,r,dvrvs,vre,dvrevs,ur,uy] = mat_disperse(thk,dns,vp,vs,freq,offsets);
toc
