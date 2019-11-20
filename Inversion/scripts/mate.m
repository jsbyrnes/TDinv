function [ o_z, o_vs, o_vpvs ] = mate( population, layers, quality, elistism, m_z, m_vs, m_vpvs )
%MATE Makes all the kiddies

    o_z = zeros(layers, population);
    o_vs = zeros(layers, population);
    o_vpvs = zeros(layers, population);

    %first sort by quality
    
    sorting = sort(quality);
            
    cutind = floor(population*elistism);
        
    %find the quality value that divides the better/worst models
    
    mark = sorting(cutind);
    
    ip = (quality <= mark);
    
    num = sum(ip);
    
    %make the top set and store their genes
    t_z = m_z(:, ip);
    t_vs = m_vs(:, ip);
    t_vpvs = m_vpvs(:, ip);
    
    %now pick cutind random ones and store their genes
    
    r_z = m_z(:, randi([1 population], 1, num));
    r_vs = m_vs(:, randi([1 population], 1, num));
    r_vpvs = m_vpvs(:, randi([1 population], 1, num));    
    
    %first dimension is layers, second dimension is member index, third
    %dimension is 1=top 2=random
    
    p_z = cat(3, t_z, r_z);
    p_vs = cat(3, t_vs, r_vs);
    p_vpvs = cat(3, t_vpvs, r_vpvs);
        
    for i = 1:num
            
        for j = 1:layers
            
            which = randi([1 2], 1, 1);
        
            o_z(j, i) = squeeze(p_z(j,i,which));
            o_vs(j, i) = squeeze(p_vs(j,i,which));
            o_vpvs(j, i) = squeeze(p_vpvs(j,i,which));

    
        end
        
    end
    
    left = population - num;
    
    r1_z = m_z(:, randi([1 population], 1, left));
    r1_vs = m_vs(:, randi([1 population], 1, left));
    r1_vpvs = m_vpvs(:, randi([1 population], 1, left));    

    r2_z = m_z(:, randi([1 population], 1, left));
    r2_vp = m_vs(:, randi([1 population], 1, left));
    r2_vpvs = m_vpvs(:, randi([1 population], 1, left));    

    p_z = cat(3, r1_z, r2_z);
    p_vs = cat(3, r1_vs, r2_vp);
    p_vpvs = cat(3, r1_vpvs, r2_vpvs);
    
    for i = 1:left
            
        for j = 1:layers
            
            which = randi([1 2],1,1);
        
            o_z(j, num + i) = squeeze(p_z(j,i,which));
            o_vs(j, num + i) = squeeze(p_vs(j,i,which));
            o_vpvs(j, num + i) = squeeze(p_vpvs(j,i,which));
    
        end
        
    end
    
    %offspring should now be all mated up

end

