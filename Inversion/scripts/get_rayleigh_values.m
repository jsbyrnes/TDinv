function [ HV, PV, GV ] = get_rayleigh_values( z, vs, vp, rho, frequency)
%This formats the model structure for disp80 and returns the values
%needed. Takes km,km/s, and g/cm^3.
    
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
    
        if min(vs)==0

            min_c = 0.1;
        
        else

            min_c = 0.1*min(vs);

        end

        max_c = 1.5*max(vs);

        [eigenvalues, ~, PV(i), GV(i), ~, ~ ] = ...
                raylei2_gateway(1/frequency(i), min_c, max_c, (max_c - min_c)/50, 2*max(z), 2*max(z), 1e-3, h, rho, vp, vs, ...
                200, 1, 1, 0, length(vp));
                           
        HV(i) = abs(eigenvalues(3,1))./abs(eigenvalues(1,1));
                
        if isnan(HV(i))

            HV(i) = 1e9;%make it not fit

        end

    end
    
    HV(HV>3)   = 3;
    PV(PV<0.1) = 0.1;
    
end

