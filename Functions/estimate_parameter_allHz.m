%get the mode of the data

clear, close all

for i = 1:9
    
    load([ './CAG4/CAG4-' num2str(i)], 'HV_mean');
    
    [mode_allbaz(i), error_all(i),~ , ~, ~] = measure_mode(HV_mean, 100);
    
end

figure(1)
errorbar(1:9, mode_allbaz, error_all);
xlabel('Hz');
ylabel('HV ratio');

HV.value     = mode_allbaz;
HV.frequency = 1:9;
HV.error     = error_all;
save('CAG4', 'HV');

