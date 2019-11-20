function [ Model ] = define_project( name )
%DEFINE_PROJECT Defines the model structure which sets the velocity model
%being used in the grid search

    switch name

        case 'simple_syn'

            %type of algorithm to use
            Model.mode = 1; %0 = grid search
                            %1 = genetic algorithm
                            
            %define base model

            %Define the starting depth of the layers in km
            Model.z = [ 0.5 1 2 ];%last depth is for half space, make large

            %in km/s
            %velocity of the incident model
            %This is vp for a P wave and vs for an SV wave. vpvs define the
            %other
            Model.vs = [ 1 1 2 ];

            %vp/vs in model in the given layer
            Model.vpvs = [ 4 3 2 ];
        
            Model.rho = nafedrake_rho(Model.vs.*Model.vpvs);
            
            %define the range that each value can vary over
            %eg is the depth is 50 and pert 30, range is 35-65
            %GAs are best when the range is large, but not so large that it
            %doesn't have time to converge

            %it's ok if the depth ranges overlap/go negative, the code will
            %adjust for that. These values will be forced at 1km
            %spacing/1km depth.

            %last must be zero
            Model.z_pert = [ 2 2 0 ];
            Model.vs_pert = [ 2 2 0 ];
            Model.vpvs_pert = [ 0 0 0 ];
                        
            Model.min_z = 0.01;
            Model.min_vs = 0.1;
            Model.min_vpvs = 2;
            
            %error_checking(Model);
            
        otherwise
            
            error('Project Name doesnt have a model associated with it');

    end
    
end