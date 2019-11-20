%%%%%%%%invert a bessel function for a PVr curve

clear, close all

name             = 'run25mday2';
starting_u       = 0.001;%starting model is uniform velocity
uneven_factor    = 1;
iterations       = 1e5;
jump_factor      = 500;
starting_f       = [ 2:2:20 ];%logspace(0,log10(15),10);
fit_thresh       = 0.07;
nmax             = 2;%test three cycles up and down
cycle_correct    = 1;

load(name);
freq_range = [1 20];
starting_u = starting_u*ones(size(starting_f));
CZ_error   = fit_thresh*ones(size(CZ));

[ u_models, fit_history ] = invert_PVr( starting_f, starting_u, ...
    uneven_factor, iterations, CZ, CZ_error, f, freq_range, arclen, jump_factor);

figure
plot(fit_history)
xlabel('Interation');
ylabel('Fit');

[~, ind]      = min(fit_history);
fit_mean      = fit_history<fit_thresh;

for i = 1:length(starting_f)

    uset = u_models(i, fit_mean);
    
    mean_model(i)    = mean(unique(uset));
    std_model(i)     = std(unique(uset));

end
    
figure
errorbar(starting_f', mean_model, std_model);
xlabel('Frequency, Hz');
ylabel('Slowness, s/m');

mean_model_cs = zeros(length(mean_model), nmax*2 + 1);
n = -nmax:nmax;

figure
hold on
for i = 1:length(n)
   
    plot(starting_f, abs(mean_model + n(i)./(arclen*starting_f)), '.')
    
    mean_model_cs(:,i) = abs(mean_model + n(i)./(arclen*starting_f));
    
    legendtext{i} = [ 'n = ' num2str(n(i)) ];
        
end
xlabel('Hz')
ylabel('Slowness')
ylim([0 0.005])
legend(legendtext)

if cycle_correct

    for i = 2:length(starting_f)

        %choose the cycle that changed the least from before
        [~, ind] = min( abs(mean_model(i-1) - mean_model_cs(i, :)) );

        mean_model(i) = mean_model_cs(i, ind);

    end
    
    figure
    errorbar(starting_f', mean_model, std_model);
    xlabel('Frequency, Hz');
    ylabel('Slowness, s/m');
    title('Cycle corrected')
    
end

PVr.frequency      = starting_f;
PVr.value_u        = mean_model;
PVr.error          = std_model;
save('PVr_20Hzfull_Correction', 'PVr');