function [ model, action, alpha, HV, PVr ] = make_newmodel_H( model, inverse_parameters, HV, PVr )
%make_newmodel Makes a new model
%1 of five things can happen. 
%1) Add a new knot to the model
%2) Remove a knot from the model
%3) Change a preexisting knot's value. Could be any of the three parameters
%on a particular knot. 
%4) Move one of the knots around
%5) Modify the uncertainties on the data (not implemented yet)

    action = 0;

    %note that inverse_parameters is not returned, so this is reset at
    %the end of the script.
%     cooldown = (1 + inverse_parameters.thermal_parameters(1)*erf(iteration/inverse_parameters.thermal_parameters(2)));
%     inverse_parameters.sigdepths = inverse_parameters.sigdepths*cooldown;
%     inverse_parameters.sigvs = inverse_parameters.sigvs*cooldown;
%     inverse_parameters.sigvpvs = inverse_parameters.sigvpvs*cooldown;
                
    modelP = model;
    HVP    = HV;
    PVrP   = PVr;
    
    %update vs or vpvs?
    mode = randsrc(1, 1, [ [1 2]; inverse_parameters.weights]);
    
    if mode == 1
       
        mode = 'vs';
        
    elseif mode == 2
        
        mode = 'vpvs';
        
    end
    
    while action == 0

        action = randi([ 1 4 ], 1);%hierchial currently turned off
        %action = 3;%for debugging
                    
        if action == 2 && (model.(mode).n == 2)%make sure there's a top and bottom at least
           
            action = 0;
            
        elseif action == 4 && (model.(mode).n == 2)%don't move the top and bottom
           
            action = 0;
            
        end
           
        if action == 1 && (inverse_parameters.max_knots <= model.(mode).n)
        
            action = 0;
            
        end
            
    end
        
    alpha = 1;
    
    switch action
        
        case 1 %birth
            
            model.(mode).n = model.(mode).n + 1;
                     
            model.(mode).z(model.(mode).n) = inverse_parameters.depth*rand(1);
                            
            %get the vs and vpvs already at that depth
            pre_value = interp1(model.(mode).z(1:(model.(mode).n - 1)), model.(mode).value(1:(model.(mode).n - 1)), ...
                model.(mode).z(model.(mode).n), inverse_parameters.interp_style, 'extrap');
            
            %make new value
            model.(mode).value(model.(mode).n) = normrnd(pre_value, inverse_parameters.sig.(mode));
            
            if model.(mode).value(model.(mode).n) < inverse_parameters.limits.(mode)(1) || model.(mode).value(model.(mode).n) > inverse_parameters.limits.(mode)(2)
            
                alpha = -1;
                
            end
                                                                            
        case 2 %death
                                    
            kill = randi(model.(mode).n, 1);%don't remove the top or bottom
                        
            model.(mode).z(kill)     = [];
            model.(mode).value(kill) = [];
            model.(mode).n           = model.(mode).n - 1;
            
            %get the new vs/vpvs at this spot
            pre_value   = interp1(model.(mode).z, model.(mode).value, ...
                modelP.(mode).z(kill), inverse_parameters.interp_style, 'extrap');
            
            if model.(mode).n < 4

                alpha = -1;

            end

        case 3 %change
            
            change = randi(model.(mode).n + 1, 1) - 1;%I want to be able to select zero
                                                
            if change == 0%static shift
                
                pert = normrnd(0, inverse_parameters.sig.(mode));
                
                model.(mode).value = model.(mode).value + pert;
                
                for kk = 1:length(model.(mode).value)
                   
                    if model.(mode).value(kk) < inverse_parameters.limits.(mode)(1) || model.(mode).value(kk) > inverse_parameters.limits.(mode)(1)

                        alpha = -1;

                    end
                                                        
                end
                
            else%individual parameter shifts
                
                model.(mode).value(change)   = normrnd(model.(mode).value(change), inverse_parameters.sig.(mode));
                
                if model.(mode).value(change) < inverse_parameters.limits.(mode)(1) || model.(mode).value(change) > inverse_parameters.limits.(mode)(2)
                    
                    alpha = -1;
                    
                end
                                
            end
            
        case 4 %move
            
            move = randi(model.(mode).n, 1);
            
            pert                 = normrnd(0, inverse_parameters.sig.depths);
            model.(mode).z(move) = model.(mode).z(move) + pert;
                        
            if model.(mode).z(move) < 0 || model.(mode).z(move) > inverse_parameters.depth

                alpha = -1;

            end

        case 5 %hyperparameters
                                                
            par = rand(1);
            
            if par < 0.5
                
                HV.error(:) = normrnd(HV.error(1), inverse_parameters.sigHV); 
                
            elseif par >= 0.5
                
                PVr.error(:) = normrnd(PVr.error(1), inverse_parameters.sigPVr); 
                
            end
                                            
    end
                       
    [ model ] = evaluate_models_matdisp(model, HV, PVr, inverse_parameters);
    
    if model.valid == 0
        
        alpha =  -1;
        
    end
    
    if alpha == 1
    
        switch action

            case 1

                alpha = min([1 (inverse_parameters.sig.(mode)*sqrt(2*pi)/diff(inverse_parameters.limits.(mode)))*...
                    exp( ((model.(mode).value(end) - pre_value)^2)/(2*inverse_parameters.sig.(mode)^2) - (model.fit - modelP.fit)/2)]);

            case 2

                if model.fit >= modelP.fit
                
                    alpha = min([1 diff(inverse_parameters.limits.(mode))/(inverse_parameters.sig.(mode)*sqrt(2*pi))*...
                        exp( -((modelP.(mode).value(kill) - pre_value)^2)/(2*inverse_parameters.sig.(mode)^2) - (model.fit - modelP.fit)/2)]);
                    
                else
                    
                    alpha = 1;
                    
                end

            case 3

                alpha = min([1 exp(-(model.fit - modelP.fit)/2) ]);

            case 4

                alpha = min([1 exp(-(model.fit - modelP.fit)/2) ]);

            case 5 %not currently being used
                
                if par < 0.5
                
                    alpha = min( [ log(1) (log((HVP.error(1)/HV.error(1)))*(length(HV.error)) - ((model.HVfit - modelP.HVfit)/2)) ] );
                
                elseif par >= 0.5
                    
                    alpha = min( [ log(1) (log((PVrP.error(1)/PVr.error(1)))*(length(PVr.error)) - ((model.PVrfit - modelP.PVrfit)/2)) ] );
    
                end
                    
                alpha = exp(alpha);%done with log domain for under/over flow errors on the lines above
                
        end
    
    end
    
    %now sort by depth
    [ model.(mode).z, ind] = sort(model.(mode).z);
    model.(mode).value     = model.(mode).value(ind);
    
end

