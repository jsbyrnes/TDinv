function [ model_fit, GVr_models ] = evaluate_models_GVr( thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m, GVr, inverse_parameters)
%evaluate_models Get GVr ratios and a fit for all of the input models

    [~, n] = size(thickness);
    
    nfreq = length(GVr.frequency);
    
    GVr_models = zeros(nfreq, n);
    
    %parfor i = 1:n
    for i = 1:n
                
        %convert the parameters to a velocity model
        [z, vs] = make_parameterized_model(thickness(:, i), vs_b(:, i), vs_m(:, i), inverse_parameters.depth, inverse_parameters.dz);
        [~, vpvs] = make_parameterized_model(thickness(:, i), vpvs_b(:, i), vpvs_m(:, i), inverse_parameters.depth, inverse_parameters.dz);
        [~, rho] = make_parameterized_model(thickness(:, i), rho_b(:, i), rho_m(:, i), inverse_parameters.depth, inverse_parameters.dz);
        
        %get the GVr ratio for the model
        [ ~, ~, GVr_models(i) ] = get_rayleigh_values(z, vs, vs.*vpvs, rho, GVr.frequency);
        
        %get the misfit relative to the error
        model_fit(i) = mean(((GVr_models(:, i) - GVr.value)./GVr.error).^2);
        
    end

end

