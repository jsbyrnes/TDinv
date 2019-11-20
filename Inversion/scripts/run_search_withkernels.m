function modelhist = run_search( HV, PVr, inverse_parameters, starting_set, chainind)
%Does the search

    disp(['Chain #' num2str(chainind) ' initiating' ]);
    rng('shuffle')
    validmodel = 0;    
    
    count = -1;
    
    while ~validmodel
    
        count = count + 1;
        
        if count > 1
            
            disp([ 'Chain #' num2str(chainind) ' restarting initialization, attempt #' num2str(count) ]);
            
        end
        
        [modeltmp] = build_startingpoints(inverse_parameters, starting_set);

        [ modeltmp.fit, modeltmp.llh, modeltmp.HVfit, modeltmp.PVrfit, modeltmp.HV, modeltmp.PVr, prior_f] ...
            = evaluate_models_matdisp(modeltmp, HV, PVr, inverse_parameters, []);
        
        if inverse_parameters.update_kernels > 0
        
            [modeltmp, validmodel, prior_f] = make_kernels(modeltmp, inverse_parameters, ...
                unique([ PVr.frequency HV.frequency]), prior_f);
            
        else
            
            modeltmp = interpolate_model(modeltmp, inverse_parameters);
            validmodel = 1;%if you aren't using kernels, this doesn't matter
            
        end
        
    end
    
    model = modeltmp;%modeltmp is a placeholder
    clear modeltmp
    
    model.fitmarks = ones(size(inverse_parameters.fit_markers));
    model.chain    = chainind;
    model.iter     = 0;
    modelhist      = model;

    modelhist = rmfield(modelhist, 'interp');

    if inverse_parameters.update_kernels > 0
       
        modelhist = rmfield(modelhist, 'kernels');
        
    end
    
    for j = 1:inverse_parameters.iter

        model.iter = j;
        
        if j == 1

            disp(['Chain #' num2str(chainind) ' has started with a fit of '...
                num2str(model.fit/(length(model.HV) + length(model.PVr))) ]);

        end

        if mod(j, inverse_parameters.writeOutOn) == 0

            disp(['Chain #' num2str(chainind) ' on iteration ' num2str(j) ...
                ' with mean fit ' num2str(model.fit/(length(model.HV) + length(model.PVr))) ]);

        end
        
        if mod(j, inverse_parameters.update_kernels) == 0 && inverse_parameters.update_kernels > 0
            
            validmodel = 0;%sometimes the model doesn't work with matdisp
            count      = 0;
            
            while ~validmodel

                [model, validmodel, prior_f] = make_kernels(model, inverse_parameters, ...
                    unique([ PVr.frequency HV.frequency]), prior_f);

                if ~validmodel
                   
                    count = count + 1;
                    
                    if length(modelhist) > count
                    
                        model = modelhist(end - count);%back up to get to valid model
                    
                        disp(['Chain #' num2str(chainind) ' backing up on iteration ' num2str(j) ]);
                        
                    else %start over :(
                        
                        model = build_startingpoints(inverse_parameters, starting_set);
                        
                        model.fitmarks = ones(size(inverse_parameters.fit_markers));
                        model.chain    = chainind;
                        model.iter     = 0;

                        [ model.fit, model.llh, model.HVfit, model.PVrfit, model.HV, model.PVr, prior_f] ...
                            = evaluate_models_matdisp(model, HV, PVr, inverse_parameters, prior_f);

                        [model, validmodel, prior_f] = make_kernels(model, inverse_parameters, ...
                            unique([ PVr.frequency HV.frequency]), prior_f);
                        
                        disp(['Chain #' num2str(chainind) ' restarting on interation ' ...
                            num2str(j) ' :(' ]);
                        
                    end
                    
                end
                
            end
            
            %update the model fit - may not be where it thinks it is...
            [ model.fit, model.llh, model.HVfit, model.PVrfit, ...
                model.HV, model.PVr] = evaluate_models_kernels( model, HV, PVr, inverse_parameters);
            
            %disp([ 'Chain #' num2str(chainind) ' updated fit from ' ...
            %    num2str(oldfit/(ndata*length(model.HV))) ' to ' num2str(model.fit/(ndata*length(model.HV))) ]);
            
        end

        alpha = -1;
        
        while alpha == -1
        
            [proposed_model, ~, alpha, prior_f] = make_newmodel(model, ...
                inverse_parameters, HV, PVr, prior_f);
            
        end

        r    =  rand(1);
        if r <= alpha

            model = proposed_model;
                        
        end
        
        if mod(j, inverse_parameters.saveEach) == 0 && j >= ...
                inverse_parameters.burnin%always save at inc. Need history regardless of "burnin"

            tmp = model;
            tmp = rmfield(tmp, 'interp');
            
            if inverse_parameters.update_kernels > 0

                tmp = rmfield(tmp, 'kernels');

            end
            
            modelhist(end+1) = tmp;
            
        end

    end
        
end
