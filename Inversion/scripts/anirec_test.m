theta = [ 0 0 0 ];
phi = [ 0 0 0 ];
z = [ 25 80 200 ];
vp = [ 6.5 8 7 ];
A = [ 0 0 0 ];
B = [ 0 0 0 ];
vpvs = [ 1.9 1.76 1.76 ];
C = [ 0 0 0 ];
rho = .33 + .77*vp;

%stairs(model.vp, model.z), title('Test Velocity Model'), grid on
%set(gca, 'YDir', 'reverse'), ylabel('Depth(km)'), xlabel('P Wave velocity(km/s)')

phase = 'SV';

baz = 0;

slow = .09;

dt = .05;

tic

[P_comp, SV_comp, ~] = anirec( phase, dt, slow, baz, vp, vp./vpvs, z, rho, A, B, C, 0 , 0, 0);

toc

% P_comp = bandpassfilt(P_comp, .05, .1, .03);
% SV_comp = bandpassfilt(SV_comp, .05, .1, .03);
% 
% [ amp, ind ] = max(SV_comp);
% 
% P_comp = P_comp/amp;
% SV_comp = SV_comp/amp;
% 
% cut = wrev(P_comp(1:ind*2));
% 
% answer = [cut;zeros(2000-length(cut), 1)];