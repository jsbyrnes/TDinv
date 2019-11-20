function [ Greens ] = make_greens_struct( Model, name)
%MAKE_GREENS_STRUCT Creates and partially fills in the descriptive part of
%the Green's Function structure
        
    %This will run if the file doesn't exist yet
    %It might take a long time make the green's functions so you want to be
    %able to just load them if you've made it before

    %Define part of the structure
    %Traces will be filled in when the code generates the green's functions
    %with the synthetics software

    %Here, count the models, then make them one by one and store them
    %in their own structure, one for each slowness

    %First, define arrays to loop over by making vectors consisting of
    %the inputs of each perturbating value

    %also does some error checking

    %z perts

    [z_parameters, z_indicies, z_pointer] = create_set(Model.z, Model.z_pert, Model.z_inc);

    %vp perts

    [vp_parameters, vp_indicies, vp_pointer] = create_set(Model.vp, Model.vp_pert, Model.vp_inc);

    %vpvs perts

    [vpvs_parameters, vpvs_indicies, vpvs_pointer] = create_set(Model.vpvs, Model.vpvs_pert, Model.vpvs_inc);

    %grad perts
    %This behaves a little differently, so its own function
    %The other parameters vary over -value to value
    %This varies from zero to value
    [grad_parameters, grad_indicies, grad_pointer] = create_set_grad(Model.grad, Model.grad_pert, Model.grad_inc);
        
    %Now fill in the Green Function    
    %There are NaNs in the parameter matricies, but if those are accessed
    %the indicies matricies are wrong anyway
        
    % I'm not entirely sure how to initialize this correctly yet
    % I suspect that would help a little with performance
    % Probably not super important
    
    %Does it need to be for loops in matlab? I think so, not sure
    %It does seems a little excessive
    
    l_rp = length(Model.rp);
    l_z = length(z_indicies);
    l_vp = length(vp_indicies);
    l_vpvs = length(vpvs_indicies);
    l_grad = length(grad_indicies);
    
    numberofmodels = l_rp*l_z*l_vp*l_vpvs*l_grad;
    
    samples = Model.len/Model.dt + 1;
    
    layers = length(Model.z);
    
    Greens(1:numberofmodels,1) = struct('z', zeros(1, layers), 'vp', zeros(1, layers), 'vpvs', zeros(1, layers), 'grad', zeros(1, layers), 'rms', [], 'rp', [], 'zerot', []);
    
    %index being used for the Greens struct
    count = 1;
    
    %so, shit, this could probably be done better
    %6 nested for loops seems... inefficient
    %However, the code doesn't take long for reasonable problem
    
    for rp_i = 1:l_rp
        
        for z_i = 1:l_z
            
            for vp_i = 1:l_vp
                
                for vpvs_i = 1:l_vpvs
                    
                    for grad_i = 1:l_grad
                                                
                        %ray parameter in s/km
                        Greens(count).rp = Model.rp(rp_i);
                        
                        %To do the more complicated values, set the model
                        %to base, then change the perturbed points only
                        
                        Greens(count).z = Model.z;
                        
                        if l_z ~= 1
                            
                            index = z_indicies(:, z_i)';
                            
                            for i = 1:length(z_pointer)
                                
                                Greens(count).z(z_pointer(i)) = z_parameters(i, index(i));
                                
                            end
                            
                        end
                        
                        %just repeat this for all 4 values getting looped
                        %Could be put in a function, but its short and very
                        %specific code probably not worth it
                        
                        Greens(count).vp = Model.vp;
                        
                        if l_vp ~= 1
                            
                            index = vp_indicies(:, vp_i)';
                            
                            for i = 1:length(vp_pointer)
                                
                                Greens(count).vp(vp_pointer(i)) = vp_parameters(i, index(i));
                                
                            end
                            
                        end
                        
                        Greens(count).vpvs = Model.vpvs;
                        
                        if l_vpvs ~= 1
                            
                            index = vpvs_indicies(:, vpvs_i)';
                            
                            for i = 1:length(vpvs_pointer)
                                
                                Greens(count).vpvs(vpvs_pointer(i)) = vpvs_parameters(i, index(i));
                                
                            end
                            
                        end
                        
                        Greens(count).grad = Model.grad;
                        
                        if l_grad ~= 1
                            
                            index = grad_indicies(:, grad_i)';
                            
                            for i = 1:length(grad_pointer)
                                
                                Greens(count).grad(grad_pointer(i)) = grad_parameters(i, index(i));
                                
                            end
                            
                        end
                            
                        count = count + 1;

                    end
                    
                end
                
            end
            
        end
        
    end
    
end