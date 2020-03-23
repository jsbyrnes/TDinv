%%
%%%%%%%%%%%%%
%This script lets you look at each of the results together
%to see if the results are reasonable before you invert
%%%%%%%%
clear
clc
close all

addpath('./Functions/')

%%
experiment_tag    = 'DiamondArray2';
radius_threshold  = 2;%m threshold for a "cluster" of radii to stack

anisotropy = 1;
nbtsp      = 50;
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

clear ZR

figure(length(Parameters.correlations) + 1)
hold on
xlabel('Frequency, Hz');
ylabel('ZR(f)');

[m,n] = size(R_mean);

for k = 1:m
        
    for kk = 1:n
        
        ZR(k).value(kk) = exp(mean(log((Z_mean{k, kk}./R_mean{k, kk}))));%exp(mean(log(Z_mean{k, kk}./R_mean{k, kk})));
        
    end
    
    ZR(k).frequency = Parameters.central_f;
    plot(ZR(k).frequency, ZR(k).value, 'o')

end

plot(Parameters.central_f, median(reshape([ZR.value], size(R_mean')), 2), 'k')

ZR(1).value = median(reshape([ZR.value], size(R_mean')), 2)';

ZR = ZR(1);

save(['./Observations/ZR-' experiment_tag ], 'ZR', '-append');

if anisotropy
   
    for k = 1:n
        
        for kk = 1:m

            zr_tmp  = [];
            azi_tmp = [];
            
            zr_tmp  = [ zr_tmp log(Z_mean{kk, k}./R_mean{kk, k}) ];%exp(mean(log(Z_mean{k, kk}./R_mean{k, kk})));
            azi_tmp = [ azi_tmp baz_hits{kk, k} ];
                            
            for nn = 1:nbtsp

                ind      = randi(length(azi_tmp), [ 1 length(azi_tmp) ]);
                fitParam = fit_AzimAniso([ azi_tmp(ind)' zr_tmp(ind)' ],'2theta',0,1);

                A0tmp(nn) = fitParam.value(1);
                A2tmp(nn) = fitParam.value(4);
                t2tmp(nn) = fitParam.value(5);

            end
        
            A0(k, kk) = mean(A0tmp);
            A2(k, kk) = mean(A2tmp);
            t2(k, kk) = mean(t2tmp);%its an angle, so this is wonky

            A0s(k, kk) = std(A0tmp);
            A2s(k, kk) = std(A2tmp);
            t2s(k, kk) = std(t2tmp);%its an angle, so this is wonky

            clear A0tmp A2tmp t2tmp
        
        end
            
    end
    
    figure(length(Parameters.correlations) + 2)
    
    for k = 1:m
    
        subplot(311)
        hold on
        errorbar(Parameters.central_f, A0(:, k), A0s(:, k))
        xlabel('Frequency, Hz');
        ylabel('ln(ZR)');
        subplot(312)
        hold on
        errorbar(Parameters.central_f, A2(:, k), A2s(:, k))
        xlabel('Frequency, Hz');
        ylabel('Amplitude, cos(2\theta)');
        subplot(313)
        hold on
        errorbar(Parameters.central_f, t2(:, k), t2s(:, k))
        xlabel('Frequency, Hz');
        ylabel('Angle, cos(2\theta)');
    
    end
    
end

save(['./Observations/ZR-' experiment_tag ], 'ZR', 'A0', 'A2', 't2', 'A0s', 'A2s', 't2s','-append');

