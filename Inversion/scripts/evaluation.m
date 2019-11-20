function [ quality ] = evaluation( members, HV)
%EVALUATION Gets the misfit for each member of the population, in terms of
%mean squared error

    [~, n] = size(members);
    
    quality = zeros(1, n);
    
    if ~iscolumn(HV.value)
        
        HV.value = HV.value';
        
    end
    
    if ~iscolumn(HV.error)
        
        HV.error = HV.error';
        
    end
    
    for i = 1:n;
        
        diff = (HV.value - members(:, i))./HV.error;
                
        quality(i) = mean(diff.^2);
        
    end

end

