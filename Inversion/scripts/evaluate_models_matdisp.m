function [ model ] = evaluate_models_matdisp( model, HV, ZJ0, inverse_parameters)
    
    model = interpolate_model(model, inverse_parameters);
    fall  = unique([ ZJ0.frequency HV.frequency]);
    
    if ~model.valid
        
        return
        
    end
    
    [PVr_model, ~, HV_model] = mat_disperse(diff(model.interp.z),model.interp.rho...
        ,model.interp.vpvs.*model.interp.vs,model.interp.vs,fall);
            
    if HV_model ~= 0
    
        model.HV  = interp1(fall, HV_model,  HV.frequency,  'pchip');
        model.PVr = interp1(fall, PVr_model, ZJ0.frequency, 'pchip');
    
        model.ZJ0 = besselj(0, (ZJ0.frequency.*2*pi*ZJ0.r)./model.PVr);
        
        model.valid = 1;
        
    else
        
        %need to be able to handle invalid models!!!!!
        disp('Invalid model')
        model.valid = 0;
        return
        
    end
        
    %get the misfit relative to the error
    
    if strcmp(inverse_parameters.H.sig_style, 'uniform')
    
        model.HVfit = sum(((model.HV - HV.value)/model.HV_error).^2);
        HVllh       = sum(-length(HV.value)*log(model.HV_error)- 0.5*model.HVfit);
       
    elseif strcmp(inverse_parameters.H.sig_style, 'fixed')
        
        model.HVfit = sum(((model.HV - HV.value)./HV.error).^2);
        HVllh       = sum(-log(HV.error)- 0.5*model.HVfit);
        
    end
    
    if any(PVr_model) ~= 0 %some models don't have a velocity, skip them. They aren't physical.
        
        if strcmp(inverse_parameters.H.sig_style, 'fixed')

            model.ZJ0fit = sum(((model.ZJ0 - ZJ0.value)./ZJ0.error).^2);
            ZJ0llh       = sum(-log(ZJ0.error)- 0.5*model.ZJ0fit);

        elseif strcmp(inverse_parameters.H.sig_style, 'uniform')

            model.ZJ0fit = sum(((model.ZJ0 - ZJ0.value)./model.ZJ0_error).^2);
            ZJ0llh       = sum(-length(ZJ0.value)*log(model.ZJ0_error)- 0.5*model.ZJ0fit);

        end
                
    else
        
        %disp('Zero PVr')
        model.ZJ0fit = 1e10;
        model.ZJ0llh = -Inf;
        
    end
    
    if model.disable_HV == 0
    
        model.fit  = model.HVfit + model.ZJ0fit;
        model.nfit = model.fit/(length(model.HV) + length(model.ZJ0));
        model.llh  = HVllh + ZJ0llh - 0.5*log(2*pi);

        if model.ZJ0fit/length(model.ZJ0) > (model.HVfit/length(model.HV))*inverse_parameters.enforce_ZJ0 ...
                && inverse_parameters.enforce_ZJ0 > 0 %check if HV is just disabled

            model.valid = 0;

        end
        
    elseif model.disable_HV == 1
    
        model.fit  = model.ZJ0fit;
        model.nfit = model.fit/(length(model.ZJ0));
        model.llh  = (ZJ0llh - 0.5*log(2*pi));
        
    end
    
end

