function [ z, value ] = make_parameterized_model( thickness, b, m, depth, dz )
%make_parameterized_model. Given some parameters, returns a model that is based on layers
%with constant gradients. b and m define a line, that is, value = b + z*m
%in the layer. Could be easily extended to more complex models. Depth
%can be arbitarily large.

    if length(b) ~= length(m)
        
        error('Check length of b and m');
        
    end
    
    if length(thickness) ~= (length(b) - 1)
        
        error('thickness should be one shorter than b and m (half space)');
        
    end
        
    thickness(end+1) = depth;
    
    z = dz:dz:depth;
    value = zeros(size(z));
    
    for i = 1:length(b)
       
        %which z are in the layer?
        if i == 1
           
            ind = find(z<=thickness(1));
            
        else
            
            ind = find((z>sum(thickness(1:i-1))) & z<=sum(thickness(1:i)));
            
        end
        
        value(ind) = b(i) + (z(ind) - z(ind(1)))*m(i);
            
    end
    
end

