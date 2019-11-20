function [ model_fit, HV_models, thickness_new,...
            vs_b_new, vpvs_b_new, rho_b_new, ...
            vs_m_new, vpvs_m_new, rho_m_new  ] = neighborhood_search(...
            thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m, inverse_parameters, HV, best_fit)
%neighborhood_search Makes a new model and set if its any good.

    for i = 1:inverse_parameters.iter

        disp([ '->Iteration ' num2str(i)]);
        [thickness_new, vs_b_new, vpvs_b_new, rho_b_new, vs_m_new, vpvs_m_new, rho_m_new] = ...
            make_newmodel(thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m, inverse_parameters);
        
        [model_fit, HV_models] = evaluate_models(thickness_new, vs_b_new, vpvs_b_new, rho_b_new, ...
            vs_m_new, vpvs_m_new, rho_m_new, HV, inverse_parameters);
        
        %recenter the neighborhood on the best one
        [~, ind] = min(model_fit);
        
        if model_fit(ind) < best_fit
        
            thickness = thickness_new(: ,ind);
            vs_b      = vs_b_new(: ,ind);
            vpvs_b    = vpvs_b_new(: ,ind);
            rho_b     = rho_b_new(: ,ind);
            vs_m      = vs_m_new(: ,ind);
            vpvs_m    = vpvs_m_new(: ,ind);
            rho_m     = rho_m_new(: ,ind);
        
            best_fit = model_fit(ind);%update the fit.
            
        end

        disp(['     Fit is ' num2str(best_fit)]);
        
    end
        
end

