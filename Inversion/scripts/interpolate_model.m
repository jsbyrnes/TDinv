function [ model ] = interpolate_model( model, inverse_parameters )

    if inverse_parameters.log_depths
    
        model.interp.z    = 2.^(linspace(0, inverse_parameters.depth, inverse_parameters.nz));
        model.vs.z        = 2.^(model.vs.z);
        model.vpvs.z      = 2.^(model.vpvs.z);
        
    else
       
        model.interp.z    = linspace(0, inverse_parameters.depth, inverse_parameters.nz);
        
    end        

    %now sort by depth - but don't save.
    [ vsz, ind ]    = sort(model.vs.z);
    vs              = model.vs.value(ind);
    [ vpvsz, ind ]  = sort(model.vpvs.z);
    vpvs            = model.vpvs.value(ind);
    
    if model.vs.n == 1
        
        model.interp.vs = model.vs.value*ones(size(model.interp.z));
        
    else
                
        model.interp.vs   = interp1(vsz, vs, model.interp.z, inverse_parameters.interp_style, 'extrap');
        
    end
    
    if model.vpvs.n == 1
        
        model.interp.vpvs = model.vpvs.value*ones(size(model.interp.z));
        
    else
        
        model.interp.vpvs   = interp1(vpvsz, vpvs, model.interp.z, inverse_parameters.interp_style, 'extrap');
    
    end
        
    %now make anything below the last now the same as the value of the last
    %node, to mimick to top of a half space
    model.interp.vs(model.interp.z > model.vs.z(end))     = model.vs.value(end);
    model.interp.vpvs(model.interp.z > model.vpvs.z(end)) = model.vpvs.value(end);
    
    if inverse_parameters.enforce_fasthalfspace
       
        if any(model.interp.vs(end) < model.interp.vs(1:end-1) ...
                | model.interp.vs(end)*model.interp.vpvs(end) ...
                < model.interp.vs(1:end-1)*model.interp.vpvs(1:end-1) )
           
            model.valid = 0;
            
        end
        
    end
    
    %%%%%%%
    %now enforce the bounds.
    
    if any( model.interp.vs < inverse_parameters.limits.vs(1)...
            | model.interp.vs > inverse_parameters.limits.vs(2)...
            | model.interp.vpvs < inverse_parameters.limits.vpvs(1)...
            | model.interp.vpvs > inverse_parameters.limits.vpvs(2) )
       
        model.valid = 0;
        
    end
        
    model.interp.rho = nafedrake_rho( model.interp.vs.*model.interp.vpvs/1000 );
            
    if inverse_parameters.log_depths
    
        model.vs.z        = log2(model.vs.z);
        model.vpvs.z      = log2(model.vpvs.z);
                
    end        
    
end


    %enforce a maximum vp and vs at the bottom of the model
%     if max(model.interp.vs.*model.interp.vpvs) ~= model.interp.vs(end)*model.interp.vpvs(end)
%         
%         model.valid = 0;
%             
%     elseif max(model.interp.vs) ~= model.interp.vs(end)
%         
%         model.valid = 0;
%                 
%     end
    
    %force a maximum curvature - don't allow things like spikes
    %use a finer grid for curvature
%     model.vs.curvature   = gradient(gradient(model.interp.vs, inverse_parameters.dz), inverse_parameters.dz);
%     model.vpvs.curvature = gradient(gradient(model.interp.vpvs, inverse_parameters.dz), inverse_parameters.dz);
%     
%     if any(abs(model.vs.curvature) > inverse_parameters.max_vs_curvature)
%        
%         model.valid = 0;
%         
%     end
%     
%     if any(abs(model.vpvs.curvature) > inverse_parameters.max_vpvs_curvature)
%        
%         model.valid = 0;
%         
%     end

