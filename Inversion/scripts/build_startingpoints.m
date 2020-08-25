function [ model ] = build_startingpoints( inverse_parameters, HV_error, PVr_error )
%build_grid Put the grid together from the gridparameters struct. Set
%single model to 1 to make the script do just one model instead of a bunch.

    %randomly select some values from within the defined range
    model.vs.n         = randi([ inverse_parameters.min_knots round(inverse_parameters.max_knots)], 1);
    model.vpvs.n       = randi([ inverse_parameters.min_knots round(inverse_parameters.max_knots)], 1);
        
    model.vs.z         = sort( inverse_parameters.depth * rand(model.vs.n, 1));
    model.vpvs.z       = sort( inverse_parameters.depth * rand(model.vpvs.n, 1));
                
    model.vs.value     = (diff(inverse_parameters.limits.vs)   * rand(model.vs.n, 1) + inverse_parameters.limits.vs(1));
    model.vpvs.value   = (diff(inverse_parameters.limits.vpvs) * rand(model.vpvs.n, 1)  + inverse_parameters.limits.vpvs(1));
    
    %put the fastest at the bottom for both vs AND vp
    [~, ind]            = max(model.vs.value);
    tmp                 = model.vs.value(end);
    model.vs.value(end) = model.vs.value(ind);
    model.vs.value(ind) = tmp;

    [~, ind]              = max(model.vpvs.value);
    tmp                   = model.vpvs.value(end);
    model.vpvs.value(end) = model.vpvs.value(ind);
    model.vpvs.value(ind) = tmp;

    if inverse_parameters.fixed_halfspace.flag

        model.vs.z(1)   = inverse_parameters.depth;
        model.vpvs.z(1) = inverse_parameters.depth;
        
        model.vs.value(1)     = inverse_parameters.fixed_halfspace.vs;
        model.vpvs.value(1)   = inverse_parameters.fixed_halfspace.vpvs;
        
    end 
    
    if strcmp(inverse_parameters.H.sig_style, 'fixed')
        
        model.HV_error  = HV_error;
        model.ZJ0_error = PVr_error;
        
    elseif strcmp(inverse_parameters.H.sig_style, 'uniform')
       
        model.HV_error  = rand(1)*(inverse_parameters.H.sigHV_limits(2) - inverse_parameters.H.sigHV_limits(1)) ...
            + inverse_parameters.H.sigHV_limits(1);
        model.ZJ0_error = rand(1)*(inverse_parameters.H.sigZJ0_limits(2) - inverse_parameters.H.sigZJ0_limits(1)) ...
            + inverse_parameters.H.sigZJ0_limits(1);
        
    elseif strcmp(inverse_parameters.H.sig_style, 'linear')
       
        model.HV_error  = rand(2, 1)*(inverse_parameters.H.sigHV_limits(2) - inverse_parameters.H.sigHV_limits(1)) ...
            + inverse_parameters.H.sigHV_limits(1);
        model.ZJ0_error = rand(2, 1)*(inverse_parameters.H.sigZJ0_limits(2) - inverse_parameters.H.sigZJ0_limits(1)) ...
            + inverse_parameters.H.sigZJ0_limits(1);
               
    end
    
    if inverse_parameters.delay_HV
        
        model.disable_HV = 1;
        
    else
        
        model.disable_HV = 0;
        
    end
    
    model.valid     = 1;
    model.converged = 0;
    
end