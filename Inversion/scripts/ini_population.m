function [ members, m_z, m_vs, m_vpvs ] = ini_population_HV( population, Model)
%ini_population Creates the initial, randomized population for the ga

    freq = [HV.frequency];

    members = zeros(length(freq), population);
    
    %should all be the same length
    layers = length(Model.z);
    
    m_z = zeros(layers, population);
    m_vs = zeros(layers, population);
    m_vpvs = zeros(layers, population);    
    
    for i = 1:population
        
        z = Model.z - Model.z_pert/2 + Model.z_pert.*rand(1, layers);
        
        %check if layers are negative depth
        
        z(z<=0) = .1;
        
        %force layers to be increasing depth
        z = sort(z);
        
        vs = Model.vs - Model.vs_pert/2 + Model.vs_pert.*rand(1, layers);
        vpvs = Model.vpvs - Model.vpvs_pert/2 + Model.vpvs_pert.*rand(1, layers);        
        
        %use the nafe drake curve for rho, take in vp
        rho = nafedrake_rho(vs.*vpvs);
        
        %get HV ratio for the current model at each frequency
        thickness = diff(Model.z)*1000;%meters
        offsets = 1000;%doesn't matter for HV ratios
        [~,~,~,~,~,~,ur,uy] = mat_disperse(thickness,rho,vs.*vpvs*1000,vs*1000,freq,offsets);

        members(:, i) = 
                
        m_z(:, i) = z;
        m_vs(:, i) = vs;
        m_vpvs(:, i) = vpvs;

    end

end