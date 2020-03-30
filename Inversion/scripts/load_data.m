function [ ZR, ZJ0 ] = load_data( inverse_parameters )
%load_data Load data. 
%Possible data are HV ratios (HV), phase velocity of rayleigh waves (PVr),
%group velocity of rayleigh waves (GVr), phase velocity of love waves
%(PVl), and group velocity of low waves (GVl)

    load([ inverse_parameters.data_path 'ZR-' inverse_parameters.data_name '.mat' ], 'ZR');
    
    load([ inverse_parameters.data_path 'SPAC-' inverse_parameters.data_name '.mat' ], 'C', 'C_error', 'Parameters', 'r');
    
    correlation = strcmp(Parameters.correlations, 'ZZ');
    [m, ~, ~] = size(C);
    
    %%%%%%make clusters of radii

    radius_vec = zeros(size(r));
    count      = 0;

    for k = 1:length(r)

        if radius_vec(k) == 0

            count                 = count + 1;        
            d                     = abs(r - r(k)) < inverse_parameters.radius_threshold;
            radius_vec(d)         = count;
            radius_cluster(count) = r(k);

        end

    end
    
    for k = 1:max(radius_vec)
        
        if sum(radius_vec==k) > 1        
        
            ZJ0(k).value     = mean(squeeze(C(radius_vec==k, correlation, :)), 1);
        
        else
           
            ZJ0(k).value     = squeeze(C(radius_vec==k, correlation, :));
            
        end
        ZJ0(k).frequency = Parameters.freq_range;
        ZJ0(k).r         = radius_cluster(k);        
        ZJ0(k).error     = linspace(0.02, 0.2, length(Parameters.freq_range));
        
    end
        
    %check dimensions
    for k = 1:length(ZR)
        
        if iscolumn(ZR(k).frequency)
            
            ZR(k).frequency = ZR(k).frequency';
            
        end
        
        if isfield(ZR(k), 'error')
        
            if iscolumn(ZR(k).error)

                ZR(k).error     = ZR(k).error';

            end
            
        else
            
            ZR(k).error = 0.15*ones(size(ZR(k).frequency));
            
        end
        
        if iscolumn(ZR(k).value)
            
            ZR(k).value     = ZR(k).value';
            
        end
        
    end
    
    for k = 1:length(ZJ0)
    
        if iscolumn(ZJ0(k).value)

            ZJ0(k).value     = ZJ0(k).value';

        end
        
        if isfield(ZJ0, 'error')
        
            if iscolumn(ZJ0(k).error)

                ZJ0(k).error     = ZJ0(k).error';

            end
                    
        end
        
        if iscolumn(ZJ0(k).frequency)

            ZJ0(k).frequency = ZJ0(k).frequency';

        end
    
    end
    
end

