function [ model, alpha ] = make_newmodel_perturb( model, inverse_parameters, HV, PVr )
%make_newmodelN Makes a new modelN
%1 of five things can happen. 
%1) Add a new knot to the modelN
%2) Remove a knot from the modelN
%3) Change a preexisting knot's value. Could be any of the three parameters
%on a particular knot. 
%4) Move one of the knots around
%5) Modify the uncertainties on the data (not implemented yet)

    action = 0;

    modelN = model;

    %update vs or vpvs?
        
    mode = randsrc(1, 1, [ [1 2]; inverse_parameters.weights]);
    
    if mode == 1
        
        mode = 'vs';
        
    elseif mode == 2
        
        mode = 'vpvs';
        
    end
                    
    while action == 0
        
        action = randi([ 1 4 ], 1);%not doing errors

        %action = 3;%for debugging
                    
        if action == 2 && (modelN.(mode).n == 2)%make sure there's a top and bottom at least
           
            action = 0;
            
        elseif action == 4 && (modelN.(mode).n == 2)%don't move the top and bottom
           
            action = 0;
            
        end
           
        if action == 1 && (inverse_parameters.max_knots <= modelN.(mode).n)
        
            action = 0;
            
        end
            
    end
        
    alpha = -1;
    
    switch action
        
        case 1 %birth
            
            modelN.(mode).n = modelN.(mode).n + 1;
                     
            modelN.(mode).z(modelN.(mode).n) = inverse_parameters.depth*rand(1);
                            
            %get the vs and vpvs already at that depth
            pre_value = interp1(modelN.(mode).z(1:(modelN.(mode).n - 1)), modelN.(mode).value(1:(modelN.(mode).n - 1)), ...
                modelN.(mode).z(modelN.(mode).n), inverse_parameters.interp_style, 'extrap');
            
            %make new value
            modelN.(mode).value(modelN.(mode).n) = normrnd(pre_value, inverse_parameters.sig.(mode));
            
            if modelN.(mode).value(modelN.(mode).n) < inverse_parameters.limits.(mode)(1) || modelN.(mode).value(modelN.(mode).n) > inverse_parameters.limits.(mode)(2)
            
                return
                
            end
            
            [ modelN.(mode).z, ind] = sort(modelN.(mode).z);
            modelN.(mode).value     = modelN.(mode).value(ind);
            
            [ modelN ] = evaluate_models_matdisp(modelN, HV, PVr, inverse_parameters);

            if ~modelN.valid
                
                return
                
            end
            
            alpha = min([1 (inverse_parameters.sig.(mode)*sqrt(2*pi)/diff(inverse_parameters.limits.(mode)))*...
                exp( ((modelN.(mode).value(end) - pre_value)^2)/(2*inverse_parameters.sig.(mode)^2) - (modelN.fit - model.fit)/2)]);
            
            r = rand(1);
            
            if r <= alpha
               
                model = modelN;
                
            end
            
        case 2 %death
                                    
            kill = randi(modelN.(mode).n, 1);%don't remove the top or bottom
                        
            modelN.(mode).z(kill)     = [];
            modelN.(mode).value(kill) = [];
            modelN.(mode).n           = modelN.(mode).n - 1;
            
            %get the new vs/vpvs at this spot
            new_value   = interp1(modelN.(mode).z, modelN.(mode).value, ...
                model.(mode).z(kill), inverse_parameters.interp_style, 'extrap');
            
            if modelN.(mode).n < 4

                return

            end

            [ modelN.(mode).z, ind] = sort(modelN.(mode).z);
            modelN.(mode).value     = modelN.(mode).value(ind);
            
            [ modelN ] = evaluate_models_matdisp(modelN, HV, PVr, inverse_parameters);

            if ~modelN.valid
                
                return
                
            end
            
            alpha = min([1 diff(inverse_parameters.limits.(mode))/(inverse_parameters.sig.(mode)*sqrt(2*pi))*...
                exp( -((model.(mode).value(kill) - new_value)^2)/(2*inverse_parameters.sig.(mode)^2) - (modelN.fit - model.fit)/2)]);
            
            r = rand(1);
            if r <= alpha
               
                model = modelN;
                
            end
            
        case 3 %change
            
            change = randi(modelN.(mode).n + 1, 1) - 1;%I want to be able to select zero
                                                
            if change == 0%static shift
                
                pert = normrnd(0, inverse_parameters.sig.(mode));
                
                modelN.(mode).value = modelN.(mode).value + pert;
                
                for kk = 1:length(modelN.(mode).value)
                   
                    if modelN.(mode).value(kk) < inverse_parameters.limits.(mode)(1) || modelN.(mode).value(kk) > inverse_parameters.limits.(mode)(1)

                        return

                    end
                                                        
                end
                
            else%individual parameter shifts
                
                modelN.(mode).value(change)   = normrnd(modelN.(mode).value(change), inverse_parameters.sig.(mode));
                
                if modelN.(mode).value(change) < inverse_parameters.limits.(mode)(1) || modelN.(mode).value(change) > inverse_parameters.limits.(mode)(2)
                    
                    return
                    
                end
                                
            end
            
            [ modelN.(mode).z, ind] = sort(modelN.(mode).z);
            modelN.(mode).value     = modelN.(mode).value(ind);
            
            [ modelN ] = evaluate_models_matdisp(modelN, HV, PVr, inverse_parameters);

            if ~modelN.valid
                
                return
                
            end
            
            alpha = min([1 exp(-(modelN.fit - model.fit)/2) ]);
            
            r = rand(1);
                 
            if r <= alpha
                
                model = modelN;
                
            end
                
        case 4 %move
            
            move = randi(modelN.(mode).n, 1);
            
            pert                  = normrnd(0, inverse_parameters.sig.depths);
            modelN.(mode).z(move) = modelN.(mode).z(move) + pert;
                        
            if modelN.(mode).z(move) < 0 || modelN.(mode).z(move) > inverse_parameters.depth

                return

            end

            [ modelN.(mode).z, ind] = sort(modelN.(mode).z);
            modelN.(mode).value     = modelN.(mode).value(ind);
            
            [ modelN ] = evaluate_models_matdisp(modelN, HV, PVr, inverse_parameters);

            if ~modelN.valid
                
                return
                
            end
            
            alpha = min([1 exp(-(modelN.fit - model.fit)/2) ]);
            
            r = rand(1);
            
            if r <= alpha
                
                model = modelN;
                
            else %delayd rejection
               
                modelN2 = model;
                
                pert                  = normrnd(0, inverse_parameters.sig.depths/5);
                modelN2.(mode).z(move) = modelN2.(mode).z(move) + pert;

                if modelN2.(mode).z(move) < 0 || modelN2.(mode).z(move) > inverse_parameters.depth

                    return

                end

                [ modelN2.(mode).z, ind] = sort(modelN2.(mode).z);
                modelN2.(mode).value     = modelN2.(mode).value(ind);
                
                [ modelN2 ] = evaluate_models_matdisp(modelN2, HV, PVr, inverse_parameters);

                if ~modelN2.valid

                    return

                end
                              
                alpha2 = min( [ 1 exp(-(modelN2.fit - modelN.fit)/2) ]);
                
                alpha = min([ 1 exp(-(modelN2.fit - model.fit)/2)*(exp( - ((5*(modelN.(mode).value(move) - modelN2.(mode).value(move)).^2)/inverse_parameters.sig.(mode))))/...
                    (exp( - (((modelN.(mode).value(move) - model.(mode).value(move)).^2)/inverse_parameters.sig.(mode))))* ...
                    ((1 - alpha2)/(1 - alpha)) ]);
                
                r = rand(1);
                
                if r <= alpha
                    
                    model = modelN2;
                    
                end
                
            end
            
        case 5 %hyperparameters
                                                
%             par = rand(1);
%             
%             if par < 0.5
%                 
%                 HV.error(:) = normrnd(HV.error(1), inverse_parameters.sigHV); 
%                 
%             elseif par >= 0.5
%                 
%                 PVr.error(:) = normrnd(Pvr.error(1), inverse_parameters.sigPVr); 
%                 
%             end
                                            
    end
                                   
end




%             if r <= alpha
%                 
%                 model = modelN;
%                 
%             else %delayed rejection
%                
%                 modelN2 = model;
%                 
%                 if change == 0%static shift
% 
%                     pert = normrnd(0, inverse_parameters.sig.(mode)/5);
% 
%                     modelN2.(mode).value = modelN2.(mode).value + pert;
% 
%                     for kk = 1:length(modelN2.(mode).value)
% 
%                         if modelN2.(mode).value(kk) < inverse_parameters.limits.(mode)(1) || modelN2.(mode).value(kk) > inverse_parameters.limits.(mode)(1)
% 
%                             return
% 
%                         end
% 
%                     end
% 
%                 else%individual parameter shifts
% 
%                     modelN2.(mode).value(change)   = normrnd(modelN2.(mode).value(change), inverse_parameters.sig.(mode)/5);
% 
%                     if modelN2.(mode).value(change) < inverse_parameters.limits.(mode)(1) || modelN2.(mode).value(change) > inverse_parameters.limits.(mode)(2)
% 
%                         return
% 
%                     end
%                     
%                 end
%                 
%                 [ modelN2.(mode).z, ind] = sort(modelN2.(mode).z);
%                 modelN2.(mode).value     = modelN2.(mode).value(ind);
%                                 
%                 [ modelN2 ] = evaluate_models_matdisp(modelN2, HV, PVr, inverse_parameters);
% 
%                 if ~modelN2.valid
% 
%                     return
% 
%                 end
%                           
%                 alpha2 = min( [ 1 exp(-(modelN2.fit - modelN.fit)/2) ]);
%                 
%                 alpha = min([ 1 exp(-(modelN2.fit - model.fit)/2)*(exp( - ((5*(modelN.(mode).value(change) - modelN2.(mode).value(change)).^2)/inverse_parameters.sig.(mode))))/...
%                     (exp( - (((modelN.(mode).value(change) - model.(mode).value(change)).^2)/inverse_parameters.sig.(mode))))* ...
%                     ((1 - alpha2)/(1 - alpha)) ]);
%                 
%                 r = rand(1);
%                 
%                 if r <= alpha
%                     
%                     model = modelN2;
%                     
%                 end
%                 
%             end
