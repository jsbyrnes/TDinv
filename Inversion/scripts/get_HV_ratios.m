function [ values ] = get_HV_ratios( z, vs, vp, rho, frequency)
%This formats the model structure for disp80 and returns the values
%needed. Takes m,m/s, and g/cm^3.

    z = z/1000;
    vs = vs/1000;
    vp = vp/1000;
    
    if ~isrow(z)
        
        z = z';
        
    end

    if ~isrow(vs)
        
        vs = vs';
        
    end

    if ~isrow(vp)
        
        vp = vp';
        
    end
    
    if ~isrow(rho)
        
        rho = rho';
        
    end
    
    if length(z) == length(vs)-1
        
        z(end + 1) = z(end)*2;
        
    end
    
    h = diff(z);
    h(end + 1) = h(end);
    
    for i = 1:length(frequency)
    
        [eigenvalues, ~, ~, ~, ~, ~ ] = ...
                raylei2_gateway(1/frequency(i), 0.9*min(vs), max(vs), (max(vs) - 0.9*min(vs))/50, max(z), max(z), 1e-3, h, rho, vp, vs, ...
                10, 1, 1, 0, length(vp));
        
        values(i) = abs(eigenvalues(3,1))./abs(eigenvalues(1,1));
    
    end
    
end

