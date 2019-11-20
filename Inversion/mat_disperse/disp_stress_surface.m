function [z,r,dr] = disp_stress_surface(freq,vr,thk,dns,cvs,cvp,Fz)

% This function calulates the displacement-stress vectors (i.e., the eigenfunctions)
% corresponding to the phase velocities (i.e., wavenumbers) contained in the vr
% matrix

% Copyright 1999 by Glenn J. Rix and Carlo G. Lai

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Establish global parameters
global NUMPOINTS LAMBDA MAXROOT

% Calculate the vector of circular frequencies
om = 2*pi*freq;

% Determine the number of layers not including the half space
N = length(thk);

% Calulate the maximum depth for determining the displacement-stress vectors
lambda_max = LAMBDA*(real(vr(:,1)).^2 + imag(vr(:,1)).^2)./real(vr(:,1))./freq;

% Initiate the depth and displacement-stress vectors and their numerical derivatives
z = zeros(NUMPOINTS,length(freq));
r = zeros(length(freq),MAXROOT,NUMPOINTS,4);
dr = zeros(length(freq),MAXROOT,NUMPOINTS,4);%jsb changed last dimension to 4

% Loop through the frequencies
for j = 1:length(freq)
   
   % Create a vector of depths
   z(:,j) = linspace(0,lambda_max(j),NUMPOINTS)';
      
   % Loop through the modes at each frequency
   index1 = find(vr(j,:));
   for m = 1:length(index1)
      
      % Calculate the wavenumber and load vector
      k = om(j)/vr(j,index1(m));
      delqz = [0 ; k*Fz/(2*pi)];
      
      % Check to see if the phase velocity is equal to the shear wave velocity
		% or compression wave velocity of one of the layers
      epsilon = 0.0001;
      while any(abs(om(j)/k-cvs)<epsilon) | any(abs(om(j)/k-cvp)<epsilon)
   		k = k * (1+epsilon);
		end   
      
      % Calculate the PSV element matrices for each layer and generalized R/T matrices
		%[e11,e12,e21,e22,du,mu,nus,nup] = psv(thk,dns,cvp,cvs,om(j),k);
		[e11,e12,e21,e22,du,mu,nus,nup] = psv_mex(thk,dns,cvp,cvs,om(j),k);
		%[td,tu,rd,ru] = modrt(e11,e12,e21,e22,du);
		[td,tu,rd,ru] = modrt_mex(e11,e12,e21,e22,du);
		[Td,Rd] = genrt(td,tu,rd,ru);
      
      % Initialize the Cd and Cu matrices
      cd = zeros(2,1,N+1);
      cu = zeros(2,1,N+1);
      
      % Calculate Cd for the first layer
      [lamd,lamu] = updown(thk,cvp,cvs,om(j),k,0,1);
      cd(:,:,1) = (e21(:,:,1) + e22(:,:,1)*lamu*Rd(:,:,1))\delqz;
      
      % Calculate Cd and Cu for the remaining layers
      for n = 1:N
         cu(:,:,n) = Rd(:,:,n)*cd(:,:,n);
         cd(:,:,n+1) = Td(:,:,n)*cd(:,:,n);
      end
      
      % Loop through the vector of depths
      %for n = 1:NUMPOINTS
      %jsb modification - don't both looping
      for n = 1:1
         
         % Determine the layer corresponding to the current depth
         index2 = find(z(n,j) <= [cumsum(thk) ; z(NUMPOINTS,j)]);
         layer = index2(1);
         
         % Calculate the up-going and down-going matrices for this depth
         [lamd,lamu] = updown(thk,cvp,cvs,om(j),k,z(n,j),layer);
         
         % Calculate the displacement-stress vector
         r(j,m,n,:) = [e11(:,:,layer) e12(:,:,layer) ; e21(:,:,layer) e22(:,:,layer)] * ...
            [lamd zeros(2) ; zeros(2) lamu] * [cd(:,:,layer) ; cu(:,:,layer)];
        
      end

   end
end
