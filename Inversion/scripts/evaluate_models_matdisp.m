function [ model ] = evaluate_models_matdisp( model, HV, ZJ0, inverse_parameters)
    
    model = interpolate_model(model, inverse_parameters);
    fall  = unique([ ZJ0.frequency HV.frequency]);
    
    if ~model.valid
        
        return
        
    end
    
    [PVr_model, ~, HV_model] = mat_disperse(diff(model.interp.z),model.interp.rho...
        ,model.interp.vpvs.*model.interp.vs,model.interp.vs,fall);
    
    if HV_model ~= 0
        
        if inverse_parameters.log_ZR

            HV_model = log(HV_model);

        end

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
    
    if ~isempty(HV.frequency)
    
        if strcmp(inverse_parameters.H.sig_style, 'uniform')

            HV_error    = model.HV_error*ones(size(HV.value));

        elseif strcmp(inverse_parameters.H.sig_style, 'linear')

            HV_error    = linspace(model.HV_error(1), model.HV_error(2), length(HV.value));

        elseif strcmp(inverse_parameters.H.sig_style, 'fixed')

            HV_error    = HV.error;

        end

        model.HVfit = sum(((model.HV - HV.value)./HV_error).^2);
       
    else
        
        model.HVfit = 0;
        
    end
    
    if ~isempty(ZJ0.frequency)
    
        if any(PVr_model) ~= 0 %some models don't have a velocity, skip them. They aren't physical.

            if strcmp(inverse_parameters.H.sig_style, 'fixed')

                ZJ0_error    = ZJ0.error;

            elseif strcmp(inverse_parameters.H.sig_style, 'linear')

                ZJ0_error    = linspace(model.ZJ0_error(1), model.ZJ0_error(2), length(ZJ0.value));

            elseif strcmp(inverse_parameters.H.sig_style, 'uniform')

                ZJ0_error    = model.ZJ0_error*ones(size(ZJ0.value));

            end

            model.ZJ0fit = sum(((model.ZJ0 - ZJ0.value)./ZJ0_error).^2);

        else

            %disp('Zero PVr')
            model.ZJ0fit = 1e10;
            model.ZJ0llh = -Inf;

        end
       
    else
        
        model.ZJ0fit = 0;
        model.ZJ0llh = 0;
        
    end
    
    if ~isempty(ZJ0.frequency) && ~isempty(HV.frequency)
    
        model.fit  = model.HVfit + model.ZJ0fit;
        model.nfit = model.fit/(length(model.HV) + length(model.ZJ0));
        model.llh  = -0.5*log(2*pi)*(length(model.HV) + length(model.ZJ0)) + ...
            2*sum(log(HV_error)) + sum(log(ZJ0_error)) - 0.5*model.fit;
        
    elseif isempty(HV.frequency)
        
        model.fit  = model.ZJ0fit;
        model.nfit = model.fit/(length(model.ZJ0));
        model.llh  = -0.5*(log(2*pi)*(length(model.ZJ0)) + ...
            sum(log(ZJ0_error)) - model.fit);
        
    elseif isempty(ZJ0.frequency)
        
        model.fit  = model.HVfit;
        model.nfit = model.fit/(length(model.HV));
        model.llh  = -0.5*(log(2*pi)*(length(model.HV)) + ...
            sum(log(HV_error)) - model.fit);
        
    end
    
end

