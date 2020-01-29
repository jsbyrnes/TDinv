%get the mode of the data

clear, close all

%load('U32A.mat', 'HV_mean', 'std_section', 'baz_array');
load('../CalTech/VTV-0.05.mat', 'HV_mean');
    
[mode_allbaz, error_all, pdf, points, ~] = measure_mode(HV_mean, 100);

figure(1)
plot(points, pdf);
xlabel('Hv ratio');
ylabel('Probability Density');
title('Probability density function at CAG4 at 8 Hz');
