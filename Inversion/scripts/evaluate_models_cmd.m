function [ model_fit, llh, HVfit, PVrfit, HV_model, PVr_model] = evaluate_models_cmd( model, HV, PVr, GVr, inverse_parameters )
%llh is log likelyhood, _cmd means the mex file is not in use

    data_counter = 0;

    model = interpolate_model(model, inverse_parameters);%not saved
    
    %theModel has km thickness, g/cm^3, km/s, km/s
    h = diff(model.interp.z);
    h(end + 1) = h(end);
    
    theModel(:,1) = h/1000;
    theModel(:,2) = model.interp.rho;
    theModel(:,3) = (model.interp.vs.*model.interp.vpvs)/1000;
    theModel(:,4) = model.interp.vs/1000;
    
    if isempty(HV)

        HVfit = 0;
        HVllh = 0;
        
    else

        data_counter = data_counter + 1;
                   
        [~,~,HV_model] = run_disper80(1./HV.frequency, theModel, 4, inverse_parameters.depth);
        
        %get the misfit relative to the error
        HVfit = mean(((HV_model - HV.value)./HV.error).^2);
        HVllh = mean(-log(HV.error)*sqrt(2*pi) - 0.5*((HV_model - HV.value)./HV.error).^2);
        
    end

    if isempty(PVr)

        PVrfit = 0;
        PVrllh = 0;
        
    else

        data_counter = data_counter + 1;
        
        [PVr_model,~,~] = run_disper80(1./PVr.frequency, theModel, 4, inverse_parameters.depth);
                        
        if isfield(PVr, 'value')%fit slowness or velocity

            PVrfit = mean(((PVr_model - PVr.value)./PVr.error).^2);
            PVrllh = mean(-log(PVr.error)*sqrt(2*pi) - 0.5*((PVr_model - PVr.value)./PVr.error).^2);
            
        elseif isfield(PVr, 'value_u')%fit slowness or velocity

            PVrfit = mean(((1./(1000*PVr_model) - PVr.value_u)./PVr.error).^2);
            PVrllh = mean(-log(PVr.error)*sqrt(2*pi) - 0.5*((1./(1000*PVr_model) - PVr.value_u)./PVr.error).^2);
            
        end
        
    end

    if isempty(GVr)

        GVrfit = 0;
        GVrllh = 0;

    else

        data_counter = data_counter + 1;

%         [model_fit_GVr, ~] = evaluate_models_GVr(thickness, vs_b, vpvs_b, rho_b, ...
%             vs_m, vpvs_m, rho_m, GVr, inverse_parameters);

    end

    model_fit = (HVfit + PVrfit + GVrfit)/data_counter;
    llh       = (HVllh + PVrllh + GVrllh)/data_counter;
        
end

