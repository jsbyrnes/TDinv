function [HV, ZJ0, model] = syn_data_linear( inverse_parameters )

    fPV = (0.5:1:19.5)';
    fHV = (0.5:1:19.5)';

    model.vs.z       = [ 0   50   100];
    model.vs.n       = length(model.vs.z);
    model.vs.value   = [ 300 1500 1500];%[ 1000 1000 1000 1000 1000  1000 1000]
    model.vpvs.z     = [ 0   50   100];
    model.vpvs.n     = length(model.vpvs.z);
    model.vpvs.value = [ 2.5 1.7  1.7];%[ 1.7 1.7 1.7 1.7 1.7  1.7  1.7];
    
    model = interpolate_model(model, inverse_parameters);
   
    [PVr.value,~,~] = mat_disperse(diff(model.interp.z),model.interp.rho...
        ,model.interp.vpvs.*model.interp.vs,model.interp.vs,fPV,1e3);

    [~,~,HV.value] = mat_disperse(diff(model.interp.z),model.interp.rho...
        ,model.interp.vpvs.*model.interp.vs,model.interp.vs,fHV,1e3);
    
    ZJ0.r         = 25;
    ZJ0.value     = besselj(0, (ZJ0.r*2*pi*fPV)./PVr.value);

    ZJ0.error                = linspace(0.02, 0.1, length(fPV))';
    ZJ0.value                = normrnd(ZJ0.value, ZJ0.error);
    ZJ0.value(ZJ0.value > 1) = 1;
    ZJ0.frequency            = fPV; 
    
    HV.value      = normrnd(HV.value, 0.075);
    HV.error      = 0.075*ones(size(HV.value));
    HV.frequency  = fHV;

end



% function [HV, PVr, model] = syn_data_buriedlayer( inverse_parameters )
% 
%     fPV = (1:1:20)';
%     fHV = (0.5:0.5:20)';
% 
%     model.vs.z       = [ 0   10  25  35  49.9 50   100];
%     model.vs.value   = [ 500 800 800 600 600  1500 1500];%[ 1000 1000 1000 1000 1000  1000 1000]
%     model.vpvs.z     = [ 0   10  25  35  49.9 50   100];
%     model.vpvs.value = [ 1.7 1.7 1.7 2.5 2.5  1.7  1.7];%[ 1.7 1.7 1.7 1.7 1.7  1.7  1.7];
%     
%     model = interpolate_model(model, inverse_parameters);
%    
%     [PVr.value,~,~] = mat_disperse(diff(model.interp.z),model.interp.rho...
%         ,model.interp.vpvs.*model.interp.vs,model.interp.vs,fPV,1e3);
%     %PVr.value = normrnd(PVr.value, 50);
%     [~,~,HV.value] = mat_disperse(diff(model.interp.z),model.interp.rho...
%         ,model.interp.vpvs.*model.interp.vs,model.interp.vs,fHV,1e3);
%     %HV.value = normrnd(HV.value, 0.1);
% %     [~,~,HV.value] = mat_disperse(diff(model.interp.z),model.interp.rho...
% %         ,model.interp.vpvs.*model.interp.vs,model.interp.vs,10,1e3);
% 
%     PVr.value = normrnd(PVr.value, 100);
%     HV.value  = normrnd(HV.value, 0.1);
% 
%     HV.error      = 0.1*ones(size(HV.value));
%     HV.frequency  = fHV;
%     PVr.error     = 100*ones(size(PVr.value));
%     PVr.frequency = fPV;
%     
% end
