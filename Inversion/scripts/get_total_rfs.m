function [ total ] = get_total_rfs( Model )
%GET_TOTAL_RFS Calculates the total number of receiver functions for a grid
%search

    total = 1;

    if any(Model.z_range)

        range_vec = ~(Model.z_range == 0);
        inc_vec = ~(Model.z_inc == 0);

        n_z = sum(Model.z_range(range_vec)/Model.z_inc(inc_vec));

        total = total*n_z;
                
    end

    if any(Model.vincident_range)

        range_vec = ~(Model.vincident_range == 0);
        inc_vec = ~(Model.vincident_inc == 0);

        n_v = sum(Model.vincident_range(range_vec)/Model.vincident_inc(inc_vec));

        total = total*n_v;
        
    end

    if any(Model.vpvs_range)

        range_vec = ~(Model.vpvs_range == 0);
        inc_vec = ~(Model.vpvs_inc == 0);

        n_vpvs = sum(Model.vpvs_range(range_vec)/Model.vpvs_inc(inc_vec));

        total = total*n_vpvs;
                
    end

end

