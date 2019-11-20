function [ rval, index, err ] = findnearest( array, value, tol)
%FINDNEAREST find the indicies of the nearest value in array to value
%within given tolerance
%   Subtracts and finds the minimum
%   Value and tol must be scalars
%   Joseph Byrnes jbyrnes@uoregon.edu

    [ m n ] = size(value);
    
    if m ~= 1 || n ~= 1
        
        error('findnearest - value should be scaler');
        
    end
    
    [ m n ] = size(tol);
    
    if ~isempty(tol)
        
        if m ~= 1 || n ~= 1
            
            error('findnearest - tol should be scaler');
            
        end
        
    end

    diff_array = abs(array - value);
    
    if isempty(tol)
    
        [err, index] = min(diff_array);
                    
    else
        
        [index tmp1] = find(diff_array < tol);
        
        err = diff_array(index);
        
%         [indicies tmp1] = find(diff_array < tol);
        
%         diff = diff_array(indicies);
        
        %now sort by diff to find where the error goes away
        
%         errmat = [diff'; indicies']';
%         
%         errmat = sortrows(errmat);
        
%         err = errmat(1, 1);
%         
%         index = errmat(1, 2);
        
    end
    
    rval = array(index);
    
end

