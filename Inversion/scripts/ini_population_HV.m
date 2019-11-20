function [ members_HV, m_z, m_vs, m_vpvs ] = ini_population_HV( population, Model, HV)
%ini_population Creates the initial, randomized population for the ga

    members_HV = zeros(length(HV.frequency), population);
    
    %should all be the same length
    layers = length(Model.z);
    
    m_z = zeros(layers, population);
    m_vs = zeros(layers, population);
    m_vpvs = zeros(layers, population);    
    
    parfor i = 1:population
        
        z = Model.z - Model.z_pert/2 + Model.z_pert.*rand(1, layers);
                
        %force layers to be increasing depth
        z = sort(z);
        
        vs = Model.vs - Model.vs_pert/2 + Model.vs_pert.*rand(1, layers);
        vpvs = Model.vpvs - Model.vpvs_pert/2 + Model.vpvs_pert.*rand(1, layers);        
        
        %check for negatives/near zeros
        z(z<=Model.min_z) = Model.min_z;
        vs(vs<=Model.min_vs) = Model.min_vs;
        vpvs(vpvs<=Model.min_vpvs) = Model.min_vpvs;

        %use the nafe drake curve for rho, take in vp
        rho = nafedrake_rho(vs.*vpvs);

        m_z(:, i) = z;
        m_vs(:, i) = vs;
        m_vpvs(:, i) = vpvs;
        
        try
        
            members_HV(:, i) = get_HV_ratios( z, vs, vpvs, rho, HV.frequency);
            
        catch
            
            members_HV(:, i) = 1000*ones(size(HV.frequency));%impossible model, make it something with a huge misfit
            
        end
            
    end

end