function [ model_fit, llh, HVfit, PVrfit, HV_model, PVr_model] = evaluate_models_kernels( model, HV, PVr, inverse_parameters)
%llh is log likelyhood
    
    model = interpolate_model(model, inverse_parameters);
    HV_model = zeros(length(HV.frequency), 1);
        
    for i = 1:length(HV.frequency)
        
        [~, find] = min(abs(model.kernels.f - HV.frequency(i)));
        
        %             HV_model(i) = model.kernels.HV0(find) + sum(model.kernels.Vs1stHV(:,i).*(model.interp.vs - model.kernels.vs0)' ...
        %                 + model.kernels.Vs2ndHV(:,i).*((model.interp.vs - model.kernels.vs0).^2)'/2 ...
        %                 + model.kernels.VpVs1stHV(:,i).*(model.interp.vpvs - model.kernels.vpvs0)' ...
        %                 + model.kernels.VpVs2ndHV(:,i).*((model.interp.vpvs - model.kernels.vpvs0).^2)'/2);
        HV_model(i) = model.kernels.HV0(find) + sum(model.kernels.Vs1stHV(:,i).*(model.interp.vs - model.kernels.vs0)' ...
            + model.kernels.VpVs1stHV(:,i).*(model.interp.vpvs - model.kernels.vpvs0)');
        
    end
    
    %get the misfit relative to the error
    HVfit = sum(((HV_model - HV.value)./HV.error).^2);
    HVllh = (sqrt(2*pi)*HV.error).^length(HV.frequency) - 0.5*HVfit;
    
    PVr_model = zeros(length(PVr.frequency), 1);
        
    for i = 1:length(PVr.frequency)
        
        [~, find] = min(abs(model.kernels.f - PVr.frequency(i)));
        
        %             PVr_model(i) = model.kernels.PVr0(find) + sum(model.kernels.Vs1stPVr(:,i).*(model.interp.vs - model.kernels.vs0)' ...
        %                 + model.kernels.Vs2ndPVr(:,i).*((model.interp.vs - model.kernels.vs0).^2)'/2 ...
        %                 + model.kernels.VpVs1stPVr(:,i).*(model.interp.vpvs - model.kernels.vpvs0)' ...
        %                 + model.kernels.VpVs2ndPVr(:,i).*((model.interp.vpvs - model.kernels.vpvs0).^2)'/2);
        PVr_model(i) = model.kernels.PVr0(find) + sum(model.kernels.Vs1stPVr(:,i).*(model.interp.vs - model.kernels.vs0)' ...
            + model.kernels.VpVs1stPVr(:,i).*(model.interp.vpvs - model.kernels.vpvs0)'); ...
            
    end
    
    PVrfit = sum(((PVr_model - PVr.value)./PVr.error).^2);
    PVrllh = (sqrt(2*pi)*PVr.error).^length(PVr.frequency) - 0.5*PVrfit;

    model_fit = (HVfit + PVrfit);
    llh       = (HVllh + PVrllh);
    
end

