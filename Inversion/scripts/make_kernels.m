function [ model, validmodel, prior_f ] = make_kernels( model, inverse_parameters, f, prior_f)
%make sensitivity kernels for Vs, VpVs, maybe rho
%these are then use to evaluate new models. f is a vector of the relevant
%frequencies you will ever look at

    model = interpolate_model(model, inverse_parameters);

    dvs   = inverse_parameters.dvs;
    dvpvs = inverse_parameters.dvpvs;

    model.kernels.f = f;
        
%     [PVr,~,HV,~,~,~,~,~, prior_f] = mat_disperse(diff(model.interp.z),model.interp.rho...
%         ,model.interp.vpvs.*model.interp.vs,model.interp.vs,f,1e3, prior_f);
         
    HV  = model.HV;
    PVr = model.PVr;

    if any(PVr==0) || length(HV) ~= length(PVr)%invalid model
       
        model = NaN;
        validmodel = 0;
        return
       
    else
        
        validmodel = 1;
        
    end
    
    for i = 1:length(model.interp.z)

        %Vs kernels
        %perturb vs
        modelP = model;
        %modelN = model;
        modelP.interp.vs(i) = model.interp.vs(i) + inverse_parameters.dvs;
        %modelN.interp.vs(i) = model.interp.vs(i) - inverse_parameters.dvs;

        [PVrP,~,HVP,~,~,~,~,~, prior_f] = mat_disperse(diff(modelP.interp.z),modelP.interp.rho...
            ,modelP.interp.vpvs.*modelP.interp.vs,modelP.interp.vs,f,1e3, prior_f);
        %[PVrN,~,HVN,~,~,~,~,~, prior_f] = mat_disperse(diff(modelN.interp.z),modelN.interp.rho...
        %    ,modelN.interp.vpvs.*modelN.interp.vs,modelN.interp.vs,f,1e3, prior_f);

        if any(PVrP==0) || length(HVP) ~= length(PVrP)%invalid model

            model = NaN;
            validmodel = 0;
            return

        else

            validmodel = 1;

        end
        
%         if any(PVrN==0) || length(HVN) ~= length(PVrN)%invalid model
% 
%             model = NaN;
%             validmodel = 0;
%             return
% 
%         else
% 
%             validmodel = 1;
% 
%         end
        
        for j = 1:length(f)

            model.kernels.Vs1stPVr(i, j) = (PVrP(j) - PVr(j))/dvs;
            %model.kernels.Vs1stPVr(i, j) = (PVrP(j) - PVrN(j))/(2*dvs);
            %model.kernels.Vs2ndPVr(i, j) = (PVrP(j) - 2*PVr(j) + PVrN(j))/(dvs*dvs);

            model.kernels.Vs1stHV(i, j) = (HVP(j) - HV(j))/dvs;
            %model.kernels.Vs1stHV(i, j) = (HVP(j) - HVN(j))/(2*dvs);
            %model.kernels.Vs2ndHV(i, j) = (HVP(j) - 2*HV(j) + HVN(j))/(dvs*dvs);

        end
        
        %VpVs kernels
        %perturb vpvs
        modelP = model;
        %modelN = model;
        modelP.interp.vpvs(i) = model.interp.vpvs(i) + dvpvs;
        %modelN.interp.vpvs(i) = model.interp.vpvs(i) - dvpvs;

        [PVrP,~,HVP,~,~,~,~,~, prior_f] = mat_disperse(diff(modelP.interp.z),modelP.interp.rho...
            ,modelP.interp.vpvs.*modelP.interp.vs,modelP.interp.vs,f,1e3, prior_f);
        %[PVrN,~,HVN,~,~,~,~,~, prior_f] = mat_disperse(diff(modelN.interp.z),modelN.interp.rho...
        %    ,modelN.interp.vpvs.*modelN.interp.vs,modelN.interp.vs,f,1e3, prior_f);
        
        if any(PVrP==0) || length(HVP) ~= length(PVrP)%invalid model

            model = NaN;
            validmodel = 0;
            return

        else

            validmodel = 1;

        end
        
%         if any(PVrN==0) || length(HVN) ~= length(PVrN)%invalid model
% 
%             model = NaN;
%             validmodel = 0;
%             return
% 
%         else
% 
%             validmodel = 1;
% 
%         end
        
        for j = 1:length(f)

            model.kernels.VpVs1stPVr(i, j) = (PVrP(j) - PVr(j))/dvpvs;
            %model.kernels.VpVs1stPVr(i, j) = (PVrP(j) - PVrN(j))/(2*dvpvs);
            %model.kernels.VpVs2ndPVr(i, j) = (PVrP(j) - 2*PVr(j) + PVrN(j))/(dvpvs*dvpvs);

            model.kernels.VpVs1stHV(i, j) = (HVP(j) - HV(j))/dvpvs;
            %model.kernels.VpVs1stHV(i, j) = (HVP(j) - HVN(j))/(2*dvpvs);
            %model.kernels.VpVs2ndHV(i, j) = (HVP(j) - 2*HV(j) + HVN(j))/(dvpvs*dvpvs);

        end
        
    end
            
    %the last one does work so repeat the value
    model.kernels.Vs1stPVr(end, :) = model.kernels.Vs1stPVr(end-1, :);
    %model.kernels.Vs2ndPVr(end, :) = model.kernels.Vs2ndPVr(end-1, :);
    model.kernels.Vs1stHV(end, :)  = model.kernels.Vs1stHV(end-1, :);
    %model.kernels.Vs2ndHV(end, :)  = model.kernels.Vs2ndHV(end-1, :);
    
    model.kernels.VpVs1stPVr(end, :) = model.kernels.VpVs1stPVr(end-1, :);
    %model.kernels.VpVs2ndPVr(end, :) = model.kernels.VpVs2ndPVr(end-1, :);
    model.kernels.VpVs1stHV(end, :)  = model.kernels.VpVs1stHV(end-1, :);
    %model.kernels.VpVs2ndHV(end, :)  = model.kernels.VpVs2ndHV(end-1, :);
    
    model.kernels.PVr0  = PVr;
    model.kernels.HV0   = HV;
    model.kernels.vs0   = model.interp.vs;
    model.kernels.vpvs0 = model.interp.vpvs;
    
end

