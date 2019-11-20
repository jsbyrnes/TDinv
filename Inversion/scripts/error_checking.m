function [ pass ] = error_checking( Model )
%ERROR_CHECKING Reviews the model structure that you've defined and
%checks for obvious errors. Makes sure that the inc evenly divides the
%perts, that the depths are ordered correcly, that the velocities are
%positive, that the gradients don't extend too far, and that the deepest
%layer is zero

%gradient over laps will be worked out later

    a = Model.z_pert./Model.z_inc;
    a(isnan(a)) = 0;
    if ~all(a == round(a))
        
        error('The z perturbations should divide evenly');
        
    end
    
    a = Model.vp_pert./Model.vp_inc;
    a(isnan(a)) = 0;
    if ~all(a == round(a))        
        error('The vp perturbations should divide evenly');
        
    end
    
    a = Model.vpvs_pert./Model.vpvs_inc;
    a(isnan(a)) = 0;
    if ~all(a == round(a))
        
        error('The vpvs perturbations should divide evenly');
        
    end
    
    a = Model.grad_pert./Model.grad_inc;
    a(isnan(a)) = 0;
    if ~all(a == round(a))
        
        error('The grad perturbations should divide evenly');
        
    end            

    for i = 1:length(Model.z) - 2 %only to this point, the last point is 0
       
        if Model.z(i) > Model.z(i+1)
            
            error('The depths should increase monotonically');
            
        end
        
    end
    
    if Model.z(end) ~= 0
        
        error('The model needs to terminate with zero, just a quirk of the synthetics software');
        
    end
    
    if Model.phase ~= 'P' || strcmpi(Model.phase, 'SV')
        
        error('Phase should be capital P or SV');
        
    end

    pass = 1;
    
end

