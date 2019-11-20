function [ members_HV, m_z, m_vs, m_vpvs ] = replacement( members_HV, quality,...
    o_z, o_vs, o_vpvs, elistism, m_z, m_vs, m_vpvs, HV, Model)
%REPLACEMENT Keep the best members, then replace the rest with offspring
        
    [~,n] = size(members_HV);
        
    %remove duplicates by setting their quality to infinity
    for i = 1:n
        
        if quality(i) > 1e9
            
            continue %already flagged
            
        end
        
        test_z = m_z(:, i);
        test_vs = m_vs(:,i);
        
        for j = 1:n
            
            if i==j
                
                continue
                
            end
            
            test1 = mean(test_z==m_z(:,j));
            test2 = mean(test_vs==m_vs(:,j));
            test = mean([test1 test2]);
            
            if test==1
                
                quality(j) = 1e10;
                
            end
            
        end
        
    end
            
    sorting = sort(quality);
    
    cutind = floor(n*elistism);
    
    %first, keep the best ones by storing them directly into the member
    %arrays. Do this by find the cut off misfit
    
    mark = sorting(cutind);
    
    ip = (quality <= mark);
    
    num = sum(ip);
        
    members_HV(:, 1:num) = members_HV(:, ip);
    m_z(:, 1:num) = m_z(:, ip);
    m_vs(:, 1:num) = m_vs(:, ip);
    m_vpvs(:, 1:num) = m_vpvs(:, ip);
    
    %now pick the best fitting offspring that are left
        
    o_z(o_z <= Model.min_z) = Model.min_z;
    o_vs(o_vs <= Model.min_vs) = Model.min_vs;
    o_vpvs(o_vpvs <= Model.min_vpvs) = Model.min_vpvs;
            
    for i = 1:n;
        
        o_z(:, i) = sort(o_z(:, i));

    end

    parfor i = 1:n
       
        %first, double check the o_z vector is ok

        rho = nafedrake_rho(o_vs(:, i).*o_vpvs(:, i));
                                    
        try
        
            offspring_HV(:, i) = get_HV_ratios(o_z(:, i), o_vs(:, i), o_vpvs(:, i), rho, HV.frequency);
            
        catch
            
            offspring_HV(:, i) = 1000*ones(size(HV.frequency));%impossible model, make it something with a huge misfit
            
        end
        
    end
        
    quality = evaluation(offspring_HV, HV);
    
    %remove duplicates by setting their quality to infinity
    for i = 1:n
        
        if quality(i) > 1e9
            
            continue %already flagged
            
        end
        
        test_z = o_z(:, i);
        test_vs = o_vs(:,i);
        
        for j = 1:n
            
            if i==j
                
                continue
                
            end
            
            test1 = mean(test_z==o_z(:,j));
            test2 = mean(test_vs==o_vs(:,j));
            test = mean([test1 test2]);
            
            if test==1
                
                quality(j) = 1e10;
                
            end
            
        end
        
    end

    
    %find the left best
    sorting = sort(quality);
        
    left = n - num;
    mark = sorting(left);
    
    ip = find((quality <= mark));
    ip = ip(1:left);

    members_HV(:, num + 1:end) = offspring_HV(:,ip);
    m_z(:, num+1:end) = o_z(:, ip);
    m_vs(:, num+1:end) = o_vs(:, ip);
    m_vpvs(:, num+1:end) = o_vpvs(:, ip);

end