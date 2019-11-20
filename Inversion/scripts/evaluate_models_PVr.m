function [ model_fit, PVr_models ] = evaluate_models_PVr( thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m, PVr, inverse_parameters)
%evaluate_models Get PVr ratios and a fit for all of the input models

    [~, n] = size(thickness);
    
    nfreq = length(PVr.frequency);
    
    PVr_models = zeros(nfreq, n);
                    
    %convert the parameters to a velocity model
    [z, vs] = make_parameterized_model(thickness, vs_b, vs_m, inverse_parameters.depth, inverse_parameters.dz);
    [~, vpvs] = make_parameterized_model(thickness, vpvs_b, vpvs_m, inverse_parameters.depth, inverse_parameters.dz);
    
    if inverse_parameters.lock_to_vp
        
        rho = nafedrake_rho(vs.*vpvs/1000);
        
    else
           
        [~, rho] = make_parameterized_model(thickness, rho_b, rho_m, inverse_parameters.depth, inverse_parameters.dz);
    
    end
        
    %get the PVr ratio for the model
    [ ~, PVr_models, ~] = get_rayleigh_values(z/1000, vs/1000, vs.*vpvs/1000, rho, PVr.frequency);
    
    %get the misfit relative to the error
        
    if isfield(PVr, 'value')%fit slowness or velocity
    
        model_fit = mean(((PVr_models - PVr.value)./PVr.error).^2);
        
    elseif isfield(PVr, 'value_u')%fit slowness or velocity
        
        model_fit = mean(((1./(1000*PVr_models) - PVr.value_u)./PVr.error).^2);
        
    end
            
end

