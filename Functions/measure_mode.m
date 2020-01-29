function [ mode_est, error, pdf, points, bw ] = measure_mode( data, varargin )
%measure_mode takes a array of data and estimates the mode by kernel
%smoothing. Bootstrapping is used to estiamte the error. Bootstrapping is
%slow. Error is the 95% confidence interval, or 2sigma.

    %hardwired min/max values. Appropriate for HV ratios. Function screws
    %up if you have a lot of data close to the edges
    min_value = 0.1;
    max_value = 7;
            
    %[pdf, points] = ksdensity(data, 0.5:0.001:4, 'censoring', censor, 'bandwidth', bw, 'kernel', 'normal', 'weights', varargin{2});
    [bw, pdf, points, ~] = kde(data, 1000, min_value, max_value);

    [~, ind_mode] = max(pdf);
    
    mode_est = points(ind_mode);

    if length(varargin) == 1
        
        bootstrap_samples = varargin{1};
        
        %points to HV
        ind_use = 1:length(data);%find(~censor);
        bootstraps = zeros(1, bootstrap_samples);
        
        if ~isempty(ind_use)
        
            %slow
            for i = 1:bootstrap_samples

                %disp(['Bootstrap sample ' num2str(i)]);

                ind_sample = randi(length(ind_use), 1, length(ind_use));

                ind_sample = ind_use(ind_sample);

                data_sample = data(ind_sample);

                if length(varargin) >= 2


                    disp('Using a function without weights right now');
                    keyboard
                    %[pdf_error, points_error] = ksdensity(data_sample, 0.5:0.001:4, 'censoring', censor, 'bandwidth', bw, 'kernel', 'normal', 'weights', varargin{2});

                else

                    [~, pdf_error, points_error, ~] = kde(data_sample, 4000, min_value, max_value);
                    %[pdf_error, points_error] = ksdensity(data_sample, 0.5:0.001:4, 'censoring', censor, 'bandwidth', bw, 'kernel', 'normal');

                end

                [~, ind_mode] = max(pdf_error);
                bootstraps(i) = points_error(ind_mode);

            end

            error = std(bootstraps);

        else

            error = [];

        end
        
    else
        
        error = [];
        
    end

end

