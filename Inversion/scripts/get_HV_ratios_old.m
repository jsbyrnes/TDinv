function [ values ] = get_HV_ratiosOFF( z, vs, vpvs, rho, frequency)
%This formats the model structure for mat_disperse and returns the values
%needed. Takes m,m/s, and g/cm^3.

    if ~isrow(z)
        
        z = z';
        
    end

    if ~isrow(vs)
        
        vs = vs';
        
    end

    if ~isrow(vpvs)
        
        vpvs = vpvs';
        
    end
    
    if ~isrow(rho)
        
        rho = rho';
        
    end
    
    if length(z) == length(vs)-1
        
        z(end + 1) = z(end)*2;
        
    end
    
    %get HV ratio for the current model at each frequency
    thickness = diff([0 z(1:end-1)]);%meters
    offsets = 1000;%doesn't matter for HV ratios
    [~,~,~,~,~,~,ur,uy] = mat_disperse(thickness,rho,vs.*vpvs,vs,frequency,offsets);

    values = abs(ur)./abs(uy);
    
end

