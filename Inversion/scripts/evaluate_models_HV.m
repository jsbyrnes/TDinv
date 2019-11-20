function [ model_fit, HV_models ] = evaluate_models_HV( thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m, HV, inverse_parameters)
%evaluate_models Get HV ratios and a fit for all of the input models

    [~, n] = size(thickness);
    
    nfreq = length(HV.frequency);
    
    HV_models = zeros(nfreq, n);
    
    %convert the parameters to a velocity model
    [z, vs] = make_parameterized_model(thickness, vs_b, vs_m, inverse_parameters.depth, inverse_parameters.dz);
    [~, vpvs] = make_parameterized_model(thickness, vpvs_b, vpvs_m, inverse_parameters.depth, inverse_parameters.dz);
    [~, rho] = make_parameterized_model(thickness, rho_b, rho_m, inverse_parameters.depth, inverse_parameters.dz);
    
    %get the HV ratio for the model
    [ HV_models, ~, ~ ] = get_rayleigh_values( z/1000, vs/1000, vs.*vpvs/1000, rho, HV.frequency);

    %get the misfit relative to the error
    model_fit = mean(((HV_models - HV.value)./HV.error).^2);

end

