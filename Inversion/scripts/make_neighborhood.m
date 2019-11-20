function [ thickness_new, vs_b_new, vpvs_b_new, rho_b_new, vs_m_new, vpvs_m_new, rho_m_new ]...
    = make_newmodel( thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m, inverse_parameters )
%make_newmodel Makes a new model

    for j = 1:inverse_parameters.layers

        if j < inverse_parameters.layers

            thickness_new(j)   = thickness(j) + normrnd(0, inverse_parameters.thickness_sigma{j}, 1, 1);

        end

        vs_b_new(j)        = vs_b(j) + normrnd(0, inverse_parameters.vs_b_sigma{j}, 1, 1);
        vpvs_b_new(j)      = vpvs_b(j) + normrnd(0, inverse_parameters.vpvs_b_sigma{j}, 1, 1);
        rho_b_new(j)       = rho_b(j) + normrnd(0, inverse_parameters.rho_b_sigma{j}, 1, 1);
        vs_m_new(j)        = vs_m(j) + normrnd(0, inverse_parameters.vs_m_sigma{j}, 1, 1);
        vpvs_m_new(j)      = vpvs_m(j) + normrnd(0, inverse_parameters.vpvs_m_sigma{j}, 1, 1);
        rho_m_new(j)       = rho_m(j) + normrnd(0, inverse_parameters.rho_m_sigma{j}, 1, 1);

    end
    
    thickness_new(thickness_new<1) = 1;
    vs_b_new(vs_b_new< 10) = 10;
    vpvs_b_new(vpvs_b_new<1) = 0.01;
    rho_b_new(rho_b_new<1) = 0.1;
        
    %redo anything that violates the monotonic condition new
    
    %vs
    
    if inverse_parameters.vs_monotonic ~= 0
                       
        for j = 1:(inverse_parameters.layers - 1)
            
            %check if the value at the base of layers is less vs at top of next
            %layer (or reversed)
            
            vs_base = vs_b_new(j) + thickness(j).*vs_m_new(j);
            vs_top  = vs_b_new(j + 1);
            
            if inverse_parameters.vs_monotonic == 1
            
                redo_ind = vs_base > vs_top;
                
            elseif inverse_parameters.vs_monotonic == -1
                
                redo_ind = vs_base < vs_top;
                
            end
            
            count = 0;
            
            while any(redo_ind)
                
                %disp([ 'redoing a vs, ' num2str(count)]);
               
                count = count + 1;
                
                if count > 500
                    
                    break
                    
                end
                
                vs_b_new(j,redo_ind)        = vs_b(j) + normrnd(0, inverse_parameters.vs_b_sigma{j}, 1, sum(redo_ind));
                vs_m_new(j,redo_ind)        = vs_m(j) + normrnd(0, inverse_parameters.vs_m_sigma{j}, 1, sum(redo_ind));
                
                vs_base = vs_b_new(j) + thickness(j).*vs_m_new(j);
                vs_top  = vs_b_new(j+1);
                
                if inverse_parameters.vs_monotonic == 1

                    redo_ind = vs_base > vs_top;

                elseif inverse_parameters.vs_monotonic == -1

                    redo_ind = vs_base < vs_top;

                end
                
            end
            
        end
        
    end

    %vpvs
    
    if inverse_parameters.vpvs_monotonic ~= 0
                       
        for j = 1:(inverse_parameters.layers - 1)
            
            %check if the value at the base of layers is less vs at top of next
            %layer (or reversed)
            
            vpvs_base = vpvs_b_new(j) + thickness(j).*vpvs_m_new(j);
            vpvs_top  = vpvs_b_new(j+1);
            
            if inverse_parameters.vpvs_monotonic == 1
            
                redo_ind = vpvs_base > vpvs_top;
                
            elseif inverse_parameters.vpvs_monotonic == -1
                
                redo_ind = vpvs_base < vpvs_top;
                
            end
            
            count = 0;
            
            while any(redo_ind)
                
                %disp([ 'redoing a vpvs, ' num2str(count)]);
               
                count = count + 1;
                
                if count > 500
                    
                    break
                    
                end
                
                vpvs_b_new(j,redo_ind)        = vpvs_b(j) + normrnd(0, inverse_parameters.vpvs_b_sigma{j}, 1, sum(redo_ind));
                vpvs_m_new(j,redo_ind)        = vpvs_m(j) + normrnd(0, inverse_parameters.vpvs_m_sigma{j}, 1, sum(redo_ind));
                
                vpvs_base = vpvs_b_new(j) + thickness(j).*vpvs_m_new(j);
                vpvs_top  = vpvs_b_new(j+1);
                
                if inverse_parameters.vpvs_monotonic == 1

                    redo_ind = vpvs_base > vpvs_top;

                elseif inverse_parameters.vpvs_monotonic == -1

                    redo_ind = vpvs_base < vpvs_top;

                end
                
            end
            
        end
        
    end
    
    %rho
    
    if inverse_parameters.rho_monotonic ~= 0
                       
        for j = 1:(inverse_parameters.layers - 1)
            
            %check if the value at the base of layers is less vs at top of next
            %layer (or reversed)
            
            rho_base = rho_b_new(j) + thickness(j).*rho_m_new(j);
            rho_top  = rho_b_new(j+1);
            
            if inverse_parameters.rho_monotonic == 1
            
                redo_ind = rho_base > rho_top;
                
            elseif inverse_parameters.rho_monotonic == -1
                
                redo_ind = rho_base < rho_top;
                
            end
            
            count = 0;
            
            while any(redo_ind)
                
                %disp([ 'redoing a rho, ' num2str(count)]);
               
                count = count + 1;
                
                if count > 500
                    
                    break
                    
                end
                
                rho_b_new(j,redo_ind)        = rho_b(j) + normrnd(0, inverse_parameters.rho_b_sigma{j}, 1, sum(redo_ind));
                rho_m_new(j,redo_ind)        = rho_m(j) + normrnd(0, inverse_parameters.rho_m_sigma{j}, 1, sum(redo_ind));
                
                rho_base = rho_b_new(j) + thickness(j).*rho_m_new(j);
                rho_top  = rho_b_new(j+1);
                
                if inverse_parameters.rho_monotonic == 1

                    redo_ind = rho_base > rho_top;

                elseif inverse_parameters.rho_monotonic == -1

                    redo_ind = rho_base < rho_top;

                end
                
            end
            
        end
        
    end
    
    %flip the sign on any gradients that have the wrong sign
    
    %vs
    if inverse_parameters.vs_monotonic == 1
    
        flip_ind = vs_m_new < 0;
        
    elseif inverse_parameters.vs_monotonic == -1
        
        flip_ind = vs_m_new > 0;
        
    end
    
    if ~isempty(flip_ind)
        
        vs_m_new(flip_ind) = vs_m_new(flip_ind)*-1;
        
    end

    %vpvs
    if inverse_parameters.vpvs_monotonic == 1
    
        flip_ind = vpvs_m_new < 0;
        
    elseif inverse_parameters.vpvs_monotonic == -1
        
        flip_ind = vpvs_m_new > 0;
        
    end
    
    if ~isempty(flip_ind)
        
        vpvs_m_new(flip_ind) = vpvs_m_new(flip_ind)*-1;
        
    end

    %rho
    if inverse_parameters.rho_monotonic == 1
    
        flip_ind = rho_m_new < 0;
        
    elseif inverse_parameters.rho_monotonic == -1
        
        flip_ind = rho_m_new > 0;
        
    end
    
    if ~isempty(flip_ind)
        
        rho_m_new(flip_ind) = rho_m_new(flip_ind)*-1;
        
    end
    
end

