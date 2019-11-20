function [ c ] = forward_coh( sta_r, sta_a, source_phi, source_amp,  f, u)

    c = zeros(size(f));
        
    for i = 1:length(source_phi)
       
        for j = 1:length(sta_r)
           
            c = c + source_amp(i)*exp(1i*2*pi*f.*u.*sta_r(j).*...
                cos(sta_a(j) - source_phi(i)));
                        
        end
        
    end

    c = c/(length(sta_r)*length(source_phi));
    
end

