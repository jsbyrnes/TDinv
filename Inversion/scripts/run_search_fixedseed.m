function [ modelhist ] = run_search_fixedseed( HV, ZJ0, inverse_parameters, chainind, starting_model)
%Does the search
%HVhist and ZJhist will be updated if and when the code becomes hierarchical

    if ~isrow(HV.frequency)
        
        HV.frequency = HV.frequency';
        HV.error     = HV.error';
        HV.value     = HV.value';

    end
    
    if ~isrow(ZJ0.frequency)
        
        ZJ0.frequency = ZJ0.frequency';
        ZJ0.error     = ZJ0.error';
        ZJ0.value     = ZJ0.value';

    end
    
    if strcmp(inverse_parameters.H.prior_distribution, 'log-uniform')
        
        inverse_parameters.H.sigHV_limits  = log(inverse_parameters.H.sigHV_limits);
        inverse_parameters.H.sigZJ0_limits = log(inverse_parameters.H.sigZJ0_limits);

    end
                
    %so that different chains don't run the exact same monte carlo simulation.....
    %rng(round(mod(chainind*now*1e12,1e3)))
    rng(chainind)

    %weights for vs verse vpvs
    inverse_parameters.weights = inverse_parameters.weights/(sum(inverse_parameters.weights));

    disp(['Chain #' num2str(chainind) ' initiating' ]);
    
    validmodel = 0;    
    
    count = -1;
    
    while ~validmodel
    
        count = count + 1;
        
        if count > 1
            
            if inverse_parameters.H.sig_style == 0
            
                disp([ 'Chain #' num2str(chainind) ' restarting initialization, attempt #' num2str(count) ]);
            
            elseif inverse_parameters.H.sig_style == 1
               
                disp([ 'Chain #' num2str(chainind) ' restarting initialization, attempt #' num2str(count) ]);
                
            end
                
        end
        
        if nargin == 4 || isempty(starting_model)
        
            model = build_startingpoints(inverse_parameters);
            
        elseif nargin == 5
            
            model = starting_model;
            
        end
                
        model = evaluate_models_matdisp(model, HV, ZJ0, inverse_parameters);
        
        validmodel = model.valid;
            
    end
        
    model.fitmarks = ones(size(inverse_parameters.fit_markers));
    model.chain    = chainind;
    model.iter     = 0;
    modelhist      = model;

    modelhist = rmfield(modelhist, 'interp'); %save the starting point
    
    for j = 1:inverse_parameters.iter

        model.iter = j;
        
        if j == 1

            if strcmp(inverse_parameters.H.sig_style, 'fixed')
            
                disp(['Chain #' num2str(chainind) ' has started with a fit of ' num2str(model.nfit, 9) ]);
            
            elseif strcmp(inverse_parameters.H.sig_style, 'uniform')
               
                disp(['Chain #' num2str(chainind) ' has started with a fit of '...
                    num2str(model.nfit) '; HV error of ' num2str(model.HV_error, 5) ...
                    ' & ZJ0 error of ' num2str(model.ZJ0_error, 5) ]);
                
            end
            
        end

        if mod(j, inverse_parameters.writeOutOn) == 0

            if strcmp(inverse_parameters.H.sig_style, 'fixed')
            
                disp(['Chain #' num2str(chainind) ' on iteration ' num2str(j) ...
                    ' with mean fit ' num2str(model.nfit) ]);
            
            elseif strcmp(inverse_parameters.H.sig_style, 'uniform')
               
                disp(['Chain #' num2str(chainind) ' on iteration ' num2str(j) ...
                    ' with mean fit ' num2str(model.nfit) '; HV error of ' num2str(model.HV_error, 5) ' & ZJ0 error of ' num2str(model.ZJ0_error, 5)]);
                
            end
            
        end
                    
        alpha = -1;
        
        while alpha == -1 %indicates a problem
        
             [modeltmp, alpha] = make_newmodel(model, inverse_parameters, HV, ZJ0);
            
        end
        
        model = modeltmp;
        
        if model.converged == 0 && model.nfit < 1
            
            model.converged = 1;
            
        end

        if mod(j, inverse_parameters.saveEach) == 0 %always save

            tmp = model;
            tmp = rmfield(tmp, 'interp');%too keep it small
                        
            modelhist(end+1) = tmp;
            %HVhist(end+1)    = HV;
            %ZJhist(end+1)   = ZJ;
            
            %hverrorhist(j)  = HV.error(1);
            %ZJerrorhist(j) = ZJ.error(1);
            
        end

        if inverse_parameters.delay_HV && model.disable_HV == 1 && model.ZJ0fit/length(model.ZJ0) < 1
    
            disp('Enabling HV data')
            model.disable_HV = 0;
            
            %now corret the fit...
            model.HVfit = sum(((model.HV - HV.value)./HV.error).^2);
            model.fit   = (model.HVfit + model.ZJ0fit);
            model.nfit  = model.fit/(length(model.HV) + length(model.ZJ0));
            
        end
        
        fit_hist(j) = model.nfit;
        nvs_hist(j) = length(model.vs.value);
        nvp_hist(j) = length(model.vpvs.value);
        
    end
       
    %modelhist(1) = [];%starting position. Save during run so that the fields match up perfectly, 
    %but should not be included in the final ensemble. 
    %HVhist(1)  = [];
    %ZJhist(1) = [];
    
end
