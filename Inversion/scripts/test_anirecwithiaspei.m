%test the anirec function on a very large model

close all

re = 6371;

model.z = 1:1:800;

R=re./exp(model.z/re);
%R = 6371-model.z;

[model.vp, model.vs, model.rho, ~, ~, ~] = iasp91(R);

r = re - model.z;

model.vp = re./r.*model.vp;
model.vs = re./r.*model.vs;

stairs(model.vp, model.z)
title('IASPEI91 Velocity Model after Earth Flattening')
ylabel('Depth(km)'), xlabel('P wave velocity(km/s)')
set(gca, 'YDir', 'reverse')
grid on
 
model.A = zeros(1, length(model.vp));
model.theta = zeros(1, length(model.vp));
model.phi = zeros(1, length(model.vp));
model.B = zeros(1, length(model.vp));
model.C = zeros(1, length(model.vp));

[p, sv, sh] = anirec('P', .05, .0684, 0, model);

t = -5:.05:95;

figure
plot(t, sv)
title('SV Greens Function'), xlabel('Time(seconds)')
figure
plot(t, p)
title('P Greens Function'), xlabel('Time(seconds)')
