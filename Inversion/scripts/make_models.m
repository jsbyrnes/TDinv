function [ model_vector ] = make_models( Model )
%make_models Takes the model vector and gets the different models you will
%test against

    %make an individual vector for each 
    
    %z perts

    [z_parameters, z_indicies, z_pointer] = create_set(Model.z, Model.z_range, Model.z_inc);

    %vp perts

    [vincident_parameters, vincident_indicies, vincident_pointer] = create_set(Model.vincident, Model.vincident_range, Model.vincident_inc);

    %vpvs perts

    [vpvs_parameters, vpvs_indicies, vpvs_pointer] = create_set(Model.vpvs, Model.vpvs_range, Model.vpvs_inc);

    %use this to increment through everything
    count = 1;
            
    l_z = length(z_indicies);
    l_vincident = length(vincident_indicies);
    l_vpvs = length(vpvs_indicies);
    
    for z_i = 1:l_z
        
        for vincident_i = 1:l_vincident
            
            for vpvs_i = 1:l_vpvs
                                
                model_vector(count).z = Model.z;
                
                if l_z ~= 1
                    
                    index = z_indicies(:, z_i)';
                    
                    for i = 1:length(z_pointer)
                        
                        model_vector(count).z(z_pointer(i)) = z_parameters(i, index(i));
                        
                    end
                    
                end
                
                %just repeat this for all 4 values getting looped
                %Could be put in a function, but its short and very
                %specific code probably not worth it
                
                model_vector(count).vincident = Model.vincident;
                
                if l_vincident ~= 1
                    
                    index = vincident_indicies(:, vincident_i)';
                    
                    for i = 1:length(vincident_pointer)
                        
                        model_vector(count).vincident(vincident_pointer(i)) = vincident_parameters(i, index(i));
                        
                    end
                    
                end
                
                model_vector(count).vpvs = Model.vpvs;
                
                if l_vpvs ~= 1
                    
                    index = vpvs_indicies(:, vpvs_i)';
                    
                    for i = 1:length(vpvs_pointer)
                        
                        model_vector(count).vpvs(vpvs_pointer(i)) = vpvs_parameters(i, index(i));
                        
                    end
                    
                end
                
                count = count + 1;
                
            end
            
        end
        
    end
                
end

