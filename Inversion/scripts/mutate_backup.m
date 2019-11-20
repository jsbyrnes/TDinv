function [ o_z, o_vs, o_vpvs ] = mutate( o_z, o_vs, o_vpvs, rate, Model )
%MUTATE Mutate the offspring. Replace rate of the elements with a random
%number

    [n, m] = size(o_z);

    nummutate = floor(m*rate);

    %make a matrix with random location and random
    %numbers for each location

    for i = 1:m
        
        locations = randi([1 n], nummutate, 1);
        mutate_z = Model.z(locations) - Model.z_pert(locations)/2 + Model.z_pert(locations).*rand(1, nummutate);
        mutate_vs = Model.vs(locations) - Model.vs_pert(locations)/2 + Model.vs_pert(locations).*rand(1, nummutate);
        mutate_vpvs = Model.vpvs(locations) - Model.vpvs_pert(locations)/2 + Model.vpvs_pert(locations).*rand(1, nummutate);

        mutate_z(mutate_z<=Model.min_z) = Model.min_z;
        mutate_vs(mutate_vs<=Model.min_vs) = Model.min_vs;
        mutate_vpvs(mutate_vpvs<=Model.min_vpvs) = Model.min_vpvs;

        
        o_z(locations, i) = mutate_z;
        o_vs(locations, i) = mutate_vs;
        o_vpvs(locations, i) = mutate_vpvs;

    end

end

