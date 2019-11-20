function [ o_z, o_vs, o_vpvs ] = mutate( o_z, o_vs, o_vpvs, rate, Model )
%MUTATE Mutate the offspring. Replace rate of the elements with a random
%number

    [n, m] = size(o_z);

    %loop over all of the parameters in all of the offspring, and randomly
    %mutate rate% of them. 
    
    for j = 1:n
    
        for i = 1:m

            mutation_occurs = rand(1,1) < rate;
                        
            if mutation_occurs

                mutate_z = o_z(j,i) - Model.z_pert(j)/2 + Model.z_pert(j).*rand(1, 1);
                mutate_vs = o_vs(j,i) - Model.vs_pert(j)/2 + Model.vs_pert(j).*rand(1, 1);
                mutate_vpvs = o_vpvs(j) - Model.vpvs_pert(j)/2 + Model.vpvs_pert(j).*rand(1, 1);

                mutate_z(mutate_z<=Model.min_z) = Model.min_z;
                mutate_vs(mutate_vs<=Model.min_vs) = Model.min_vs;
                mutate_vpvs(mutate_vpvs<=Model.min_vpvs) = Model.min_vpvs;

                o_z(j,i) = mutate_z;
                o_vs(j, i) = mutate_vs;
                o_vpvs(j, i) = mutate_vpvs;

            end

        end

    end
end

