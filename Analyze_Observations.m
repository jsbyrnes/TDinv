%%
%%%%%%%%%%%%%
%This script lets you look at each of the results together
%to see if the results are reasonable before you invert
%%%%%%%%
clear
clc
close all
%%
experiment_tag    = 'DiamondArray';
radius_threshold  = 2;%m threshold for a "cluster" of radii to stack

%%
load(['./Observations/SPAC-' experiment_tag ]);

%%%%%%make clusters of radii

radius_vec = zeros(size(r));
count      = 0;

for k = 1:length(r)
   
    if radius_vec(k) == 0
       
        count             = count + 1;        
        d                 = abs(r - r(k)) < radius_threshold;
        radius_vec(d)     = count;
        radius_cluster(count) = r(k);
        
    end
    
end

n_radii = max(radius_vec);

[m,n,p] = size(C);

for correlation = 1:n
    
    figure(correlation)

    for k = 1:max(radius_vec)
    
        subplot(max(radius_vec),1,k)
        hold on
        
        plot(Parameters.freq_range, squeeze(C(radius_vec==k, correlation, :)), 'k--');
        plot(Parameters.freq_range, mean(squeeze(C(radius_vec==k, correlation, :)), 1), 'k', 'LineWidth', 2);        
    
        xlabel('Frequency, Hz');
        title([Parameters.correlations{correlation} ' at ' num2str(radius_cluster(k)) ' m ' ]);
        
    end

end
    
%%
%%%%%Now do ZR

load(['./Observations/ZR-' experiment_tag ]);

figure(length(Parameters.correlations) + 1)
hold on
xlabel('Frequency, Hz');
ylabel('ZR(f)');

[m,n] = size(R_mean);

for k = 1:m
    
    tmp = [];
    
    for kk = 1:n
        
        ZR(k).value(kk) = exp(mean(log((Z_mean{k, kk}./R_mean{k, kk}))));%exp(mean(log(Z_mean{k, kk}./R_mean{k, kk})));
        
    end
    
    ZR(k).frequency = Parameters.central_f;
    plot(ZR(k).frequency, ZR(k).value, 'o')

end

plot(Parameters.central_f, mean(reshape([ZR.value], size(R_mean')), 2), 'k')

save(['./Observations/ZR-' experiment_tag ], 'ZR', '-append');
