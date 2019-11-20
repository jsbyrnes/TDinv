function [ modelhist ] = run_search( HV, PVr, inverse_parameters, chainind, starting_model)
%Does the search
%HVhist and PVrhist will be updated if and when the code becomes hierarchical

    if ~isrow(HV.frequency)
        
        HV.frequency = HV.frequency';
        HV.error     = HV.error';
        HV.value     = HV.value';

    end
    
    if ~isrow(PVr.frequency)
        
        PVr.frequency = PVr.frequency';
        PVr.error     = PVr.error';
        PVr.value     = PVr.value';

    end
    
    %so that different chains don't run the exact same monte carlo simulation.....
    rng(round(mod(now*1e12,1e3)))

    %weights for vs verse vpvs
    inverse_parameters.weights = inverse_parameters.weights/(sum(inverse_parameters.weights));

    disp(['Chain #' num2str(chainind) ' initiating' ]);
    
    validmodel = 0;    
    
    count = -1;
    
    while ~validmodel
    
        count = count + 1;
        
        if count > 1
            
            disp([ 'Chain #' num2str(chainind) ' restarting initialization, attempt #' num2str(count) ]);
            
        end
        
        if nargin == 4 || isempty(starting_model)
        
            model = build_startingpoints(inverse_parameters);
        
        elseif nargin == 5
            
            model = starting_model;
            
        end
        
        model = evaluate_models_matdisp(model, HV, PVr, inverse_parameters);
        
        validmodel = model.valid;
            
    end
        
    model.fitmarks = ones(size(inverse_parameters.fit_markers));
    model.chain    = chainind;
    model.iter     = 0;
    modelhist      = model;
    %HVhist         = HV;
    %PVrhist        = PVr;

    modelhist = rmfield(modelhist, 'interp');
    
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
                    
        alpha = -1;
        
        while alpha == -1
        
             [proposed_model, ~, alpha] = make_newmodel(model, ...
                 inverse_parameters, HV, PVr);
%            [proposed_model, ~, alpha, proposed_HV, proposed_PVr] = make_newmodel_H(model, ...
%                inverse_parameters, HV, PVr);
            
        end

        r    =  rand(1);
        if r <= alpha

            model = proposed_model;
            %HV    = proposed_HV;
            %PVr   = proposed_PVr;
            
        end
        
        if mod(j, inverse_parameters.saveEach) == 0 && j > inverse_parameters.burnin

            tmp = model;
            tmp = rmfield(tmp, 'interp');%too keep it small
                        
            modelhist(end+1) = tmp;
            %HVhist(end+1)    = HV;
            %PVrhist(end+1)   = PVr;
            
            %hverrorhist(j)  = HV.error(1);
            %pvrerrorhist(j) = PVr.error(1);
            
        end

        fit_hist(j) = model.fit/(length(model.HV) + length(model.PVr));
        
    end
       
    modelhist(1) = [];%starting position. Save during run so that the fields match up perfectly, 
    %but should not be included in the final ensemble. 
    %HVhist(1)  = [];
    %PVrhist(1) = [];
    
end
