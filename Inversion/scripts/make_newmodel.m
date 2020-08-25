function [ model, alpha ] = make_newmodel( model, inverse_parameters, HV, ZJ0 )
%make_newmodelN Makes a new modelN
%1 of five things can happen. 
%1) Add a new knot to the modelN
%2) Remove a knot from the modelN
%3) Change a preexisting knot's value. Could be any of the three parameters
%on a particular knot. 
%4) Move one of the knots around
%5) Modify the uncertainties on the data (not implemented yet)

    max_count = 10;

    action = 0;
        
    %update vs or vpvs?
        
    %mode = randsrc(1, 1, [ [1 2]; inverse_parameters.weights]);
    mode = randi(2, 1);
    
    if mode == 1
        
        mode = 'vs';
        
    elseif mode == 2
        
        mode = 'vpvs';
        
    end
                    
    while action == 0
        
        if strcmp(inverse_parameters.H.sig_style, 'fixed')
        
            action = randi([ 1 4 ], 1);
                            
        else %inverting for errors
           
            action = randi([ 1 5 ], 1);%not enabled
            
        end
                    
        if action == 2 && (model.(mode).n == inverse_parameters.min_knots)
           
            action = 0;
            
        elseif action == 4 && (model.(mode).n == 1)
           
            action = 0;
            
        elseif action == 5 && (model.(mode).n == 1)
           
            action = 0;
                        
        end
           
        if action == 1 && (inverse_parameters.max_knots == model.(mode).n)
        
            action = 0;
            
        end
            
    end
        
    alpha         = -1;
    modelN        = model;%dummy step!
    modelN.valid  = 0;

    count = 0;
            
    if action ~= 4
        
        sig = inverse_parameters.sig*range(inverse_parameters.limits.(mode))/100;
        
    elseif action == 4
        
        sig = inverse_parameters.sig.*inverse_parameters.depth/100;
        
    end
            
    switch action
        
        case 1 %birth
                        
            while ~modelN.valid %techincally could be infinite
            
                modelN       = model; %will loop, and don't want to build up pert
                modelN.valid = 1;
                
                modelN.(mode).n = modelN.(mode).n + 1;
                modelN.(mode).z(modelN.(mode).n) = inverse_parameters.depth*rand(1);

                %make new value
                if strcmp(inverse_parameters.birth_style, 'perturb')

                    pre_value = interp1(model.interp.z, model.interp.(mode), ...
                        modelN.(mode).z(modelN.(mode).n), inverse_parameters.interp_style, 'extrap');

                    modelN.(mode).value(modelN.(mode).n) = normrnd(pre_value, sig);

                    modelN = interpolate_model(modelN, inverse_parameters);
                    
                    new_value = interp1(modelN.interp.z, modelN.interp.(mode), ...
                        modelN.(mode).z(modelN.(mode).n), inverse_parameters.interp_style, 'extrap');
                    
                elseif strcmp(inverse_parameters.birth_style, 'prior')

                    modelN.(mode).value(modelN.(mode).n) = range(inverse_parameters.limits.(mode))*rand(1);

                end

                if modelN.(mode).value(modelN.(mode).n) < inverse_parameters.limits.(mode)(1) || modelN.(mode).value(modelN.(mode).n) > inverse_parameters.limits.(mode)(2)

                    modelN.valid = 0;

                end

                %[ modelN.(mode).z, ind] = sort(modelN.(mode).z);
                %modelN.(mode).value     = modelN.(mode).value(ind);

                if modelN.valid
                
                    [ modelN ] = evaluate_models_matdisp(modelN, HV, ZJ0, inverse_parameters);

                end
                
                count = count + 1;
                
                if count > max_count
                    
                    return
                    
                end
                
            end
            
            if strcmp(inverse_parameters.birth_style, 'perturb')%not set up at the moment
            
                alpha = min([1 (sig*sqrt(2*pi)/diff(inverse_parameters.limits.(mode)))*...
                    exp( ((new_value - pre_value)^2)/(2*sig^2) - (modelN.fit - model.fit)/2)]);
                %alpha = min([1 (sig*sqrt(2*pi)/diff(inverse_parameters.limits.(mode)))*...
                %    exp( ((old_value - pre_value)^2)/(2*sig^2) - (modelN.fit - model.fit)/2)]);
            
            elseif strcmp(inverse_parameters.birth_style, 'prior')
               
                alpha = min([1 ((model.(mode).n + 1)/(model.(mode).n + 2))*exp(-(modelN.fit - model.fit)/2)]);
                %alpha = min([1 ((model.(mode).n)/(model.(mode).n + 2))*exp(-(modelN.fit - model.fit)/2)]);%log uniform
                
            end
                
            r = rand(1);
            
            if r <= alpha
               
                model = modelN;
                                
            end
            
        case 2 %death

            while ~modelN.valid %techincally could be infinite
            
                modelN       = model; %will loop, and don't want to build up pert
                modelN.valid = 1;
                
                if ~inverse_parameters.fixed_halfspace.flag
                
                    kill                      = randi(modelN.(mode).n, 1);
                    
                else
                   
                    kill                      = randi([ 2 modelN.(mode).n ], 1);
                    
                end
                    
                modelN.(mode).z(kill)     = [];
                modelN.(mode).value(kill) = [];
                modelN.(mode).n           = modelN.(mode).n - 1;

                %get the new vs/vpvs at this spot
                
                modelN = interpolate_model(modelN, inverse_parameters);
                
                new_value = interp1(modelN.interp.z, modelN.interp.(mode), ...
                    model.(mode).z(kill), inverse_parameters.interp_style, 'extrap');
                old_value = interp1(model.interp.z, model.interp.(mode), ...
                    model.(mode).z(kill), inverse_parameters.interp_style, 'extrap');
                
                %new_value   = interp1(modelN.(mode).z, modelN.(mode).value, ...
                %    model.(mode).z(kill), inverse_parameters.interp_style, 'extrap');
                
                if modelN.(mode).n < inverse_parameters.min_knots

                    modelN.valid = 0;

                end

                %[ modelN.(mode).z, ind] = sort(modelN.(mode).z);
                %modelN.(mode).value     = modelN.(mode).value(ind);

                if modelN.valid
                
                    [ modelN ] = evaluate_models_matdisp(modelN, HV, ZJ0, inverse_parameters);

                end
                
                count = count + 1;
                
                if count > max_count
                    
                    return
                    
                end
                    
            end
            
            if strcmp(inverse_parameters.birth_style, 'perturb')
                            
                alpha = min([1 diff(inverse_parameters.limits.(mode))/(sig*sqrt(2*pi))*...
                    exp( -((old_value - new_value)^2)/(2*sig^2) - (modelN.fit - model.fit)/2)]);
                %alpha = min([1 diff(inverse_parameters.limits.(mode))/(sig*sqrt(2*pi))*...
                %    exp( -((old_value - new_value)^2)/(2*sig^2) - (modelN.fit - model.fit)/2)]);
            
            elseif strcmp(inverse_parameters.birth_style, 'prior')
               
                alpha = min([1 (((model.(mode).n + 1)/(model.(mode).n)))*exp(-(modelN.fit - model.fit)/2)]);
                %alpha = min([1 (((model.(mode).n + 1)/(model.(mode).n))^2)*exp(-(modelN.fit - model.fit)/2)]);%log-uniform
                
            end
            
            r = rand(1);
            if r <= alpha
               
                model = modelN;
                                
            end
            
        case 3 %change
            
            while ~modelN.valid %techincally could be infinite
            
                modelN       = model; %will loop, and don't want to build up pert
                modelN.valid = 1;

                if ~inverse_parameters.fixed_halfspace.flag

                    change = randi(modelN.(mode).n, 1);
                    
                else
                   
                    change = randi([ 2 modelN.(mode).n ], 1);
                    
                end
                
                modelN.(mode).value(change)   = normrnd(modelN.(mode).value(change), sig);
                
                if modelN.(mode).value(change) < inverse_parameters.limits.(mode)(1) || modelN.(mode).value(change) > inverse_parameters.limits.(mode)(2)
                    
                    modelN.valid = 0;
                    
                end
                
                %[ modelN.(mode).z, ind] = sort(modelN.(mode).z);
                %modelN.(mode).value     = modelN.(mode).value(ind);

                if modelN.valid
                
                    [ modelN ] = evaluate_models_matdisp(modelN, HV, ZJ0, inverse_parameters);

                end
                    
                count = count + 1;
                
                if count > max_count
                    
                    return
                    
                end
                
            end
                
            alpha = min([1 exp(-(modelN.fit - model.fit)/2) ]);
            
            r = rand(1);
            
            if r <= alpha
                
                model = modelN;
                
            else %delayed rejection
                
                modelN2       = modelN; %will loop, and don't want to build up pert
                modelN2.valid = 0;
                
                while ~modelN2.valid
                
                    modelN2       = modelN; %will loop, and don't want to build up pert
                    modelN2.valid = 1;

                    if ~inverse_parameters.fixed_halfspace.flag

                        change2 = randi(modelN2.(mode).n, 1);

                    else

                        change2 = randi([ 2 modelN2.(mode).n ], 1);

                    end
                    
                    modelN2.(mode).value(change2)   = normrnd(modelN.(mode).value(change2), sig/5);

                    if modelN2.(mode).value(change2) < inverse_parameters.limits.(mode)(1) || modelN2.(mode).value(change2) > inverse_parameters.limits.(mode)(2)

                        modelN2.valid = 0;

                    end

                    %[ modelN2.(mode).z, ind] = sort(modelN2.(mode).z);
                    %modelN2.(mode).value     = modelN2.(mode).value(ind);

                    if modelN2.valid

                        [ modelN2 ] = evaluate_models_matdisp(modelN2, HV, ZJ0, inverse_parameters);

                    end

                    count = count + 1;

                    if count > max_count

                        return

                    end
                
                end
                
                q1 = (1/(sig*sqrt(2*pi)))*exp( -((modelN.(mode).value(change) - model.(mode).value(change))^2)/(2*sig^2));
                q2 = (1/(0.2*sig*sqrt(2*pi)))*exp( -((modelN2.(mode).value(change2) - modelN.(mode).value(change2))^2)/(2*(0.2*sig)^2));
                
                a1 = (1 - min([ 1 exp(-(modelN.fit - model.fit)/2) ]));
                a2 = (1 - min([ 1 exp(-(modelN.fit - modelN2.fit)/2) ]));
                
                alpha2 = min([1 exp(-(modelN2.fit - model.fit)/2)*(q2/q1)*(a2/a1) ]);%A5 in Bodin and Sambridge 2009
                          
                r = rand(1);

                if r <= alpha2

                    model = modelN2;

                end
                
            end
                                
        case 4 %move
            
            while ~modelN.valid %techincally could be infinite
            
                modelN       = model; %will loop, and don't want to build up pert
                modelN.valid = 1;
                
                if ~inverse_parameters.fixed_halfspace.flag

                    move = randi(modelN.(mode).n, 1);
                    
                else
                   
                    move = randi([ 2 modelN.(mode).n ], 1);
                    
                end
                
                pert                  = normrnd(0, sig);                
                modelN.(mode).z(move) = modelN.(mode).z(move) + pert;

                if modelN.(mode).z(move) < 0 || modelN.(mode).z(move) > inverse_parameters.depth

                    modelN.valid = 0;

                end

                %[ modelN.(mode).z, ind] = sort(modelN.(mode).z);
                %modelN.(mode).value     = modelN.(mode).value(ind);

                if modelN.valid
                
                    [ modelN ] = evaluate_models_matdisp(modelN, HV, ZJ0, inverse_parameters);

                end
                
                count = count + 1;
                
                if count > max_count
                    
                    return
                    
                end
                
            end
                            
            alpha = min([1 exp(-(modelN.fit - model.fit)/2) ]);
            
            r = rand(1);
            
            if r <= alpha
                
                model = modelN;
                
            else %delayed rejection
                
                modelN2       = modelN; %will loop, and don't want to build up pert
                modelN2.valid = 0;
                
                while ~modelN2.valid
                
                    modelN2       = modelN; %will loop, and don't want to build up pert
                    modelN2.valid = 1;

                    if ~inverse_parameters.fixed_halfspace.flag

                        move2 = randi(modelN2.(mode).n, 1);

                    else

                        move2 = randi([ 2 modelN2.(mode).n ], 1);

                    end

                    modelN2.(mode).z(move2)   = normrnd(modelN.(mode).z(move2), sig/5);

                    if modelN2.(mode).z(move2) < 0 || modelN2.(mode).z(move2) > inverse_parameters.depth

                        modelN2.valid = 0;

                    end

                    %[ modelN2.(mode).z, ind] = sort(modelN2.(mode).z);
                    %modelN2.(mode).value     = modelN2.(mode).value(ind);

                    if modelN2.valid

                        [ modelN2 ] = evaluate_models_matdisp(modelN2, HV, ZJ0, inverse_parameters);

                    end

                    count = count + 1;

                    if count > max_count

                        return

                    end
                
                end
                
                q1 = (1/(sig*sqrt(2*pi)))*exp( -((modelN.(mode).z(move) - model.(mode).z(move))^2)/(2*sig^2));
                q2 = (1/(0.2*sig*sqrt(2*pi)))*exp( -((modelN2.(mode).z(move2) - modelN.(mode).z(move2))^2)/(2*(0.2*sig)^2));
                
                a1 = (1 - min([ 1 exp(-(modelN.fit - model.fit)/2) ]));
                a2 = (1 - min([ 1 exp(-(modelN.fit - modelN2.fit)/2) ]));
                
                alpha2 = min([1 exp(-(modelN2.fit - model.fit)/2)*(q2/q1)*(a2/a1) ]);%A5 in Bodin and Sambridge 2009
                          
                r = rand(1);

                if r <= alpha2

                    model = modelN2;

                end
                
            end
                   
        case 5
            
            r = log(rand(1));
            
            if isempty(ZJ0.frequency)

                choose = 1;

            elseif isempty(HV.frequency)

                choose = 2;

            else

                choose = randi(2, 1);

            end

            if strcmp(inverse_parameters.H.sig_style, 'invert')

                HV_error    = model.HV_error*ones(size(HV.value));
                ZJ0_error   = linspace(model.ZJ0_error(1), model.ZJ0_error(2), length(ZJ0.value));
                                
                if choose == 1
                   
                    modelN.HV_error = normrnd(modelN.HV_error, inverse_parameters.H.sigHV);

                    if modelN.HV_error < inverse_parameters.H.sigHV_limits(1) || modelN.HV_error > inverse_parameters.H.sigHV_limits(2)
                       
                        r = Inf;
                        
                    end
                    
                elseif choose == 2
                    
                    choose = randi(2, 1);
                    modelN.ZJ0_error(choose) = normrnd(modelN.ZJ0_error(choose), inverse_parameters.H.sigZJ0);
                    
                    if modelN.ZJ0_error(choose) < inverse_parameters.H.sigZJ0_limits(1) || modelN.ZJ0_error(choose) > inverse_parameters.H.sigZJ0_limits(2)
                       
                        r = Inf;
                        
                    end
                    
                end
                
                HV_errorN   = modelN.HV_error*ones(size(HV.value));
                ZJ0_errorN  = linspace(model.ZJ0_error(1), model.ZJ0_error(2), length(ZJ0.value));

            end
            
            if ~isempty(ZJ0.frequency) && ~isempty(HV.frequency)
            
                %get fits without running mat_disp, makes this a free step
                modelN.ZJ0fit = sum(((modelN.ZJ0 - ZJ0.value)./ZJ0_errorN).^2);
                modelN.HVfit  = sum(((modelN.HV - HV.value)./HV_errorN).^2);
                modelN.fit    = modelN.ZJ0fit + modelN.HVfit;
                modelN.nfit   = modelN.fit/(length(HV.value) + length(ZJ0.value));
                modelN.llh    = -0.5*(log(2*pi)*(length(modelN.HV) + length(modelN.ZJ0)) + ...
                    sum(log(HV_errorN)) + sum(log(ZJ0_errorN)) - modelN.fit);
                %uncorrelated errors, so the determinent is just the product of
                %the trace, here as a log
                detold = sum(log(HV_error))  + sum(log(ZJ0_error));
                detnew = sum(log(HV_errorN)) + sum(log(ZJ0_errorN));
            
            elseif isempty(ZJ0.frequency)
                
                %get fits without running mat_disp, makes this a free step
                modelN.HVfit  = sum(((modelN.HV - HV.value)./HV_errorN).^2);
                modelN.fit    = modelN.HVfit;
                modelN.nfit   = modelN.fit/(length(HV.value));
                modelN.llh    = -0.5*(log(2*pi)*(length(modelN.HV)) + ...
                    sum(log(HV_errorN)) - modelN.fit);
                %uncorrelated errors, so the determinent is just the product of
                %the trace, here as a log
                detold = sum(log(HV_error));
                detnew = sum(log(HV_errorN));
                
            elseif isempty(HV.frequency)
                
                %get fits without running mat_disp, makes this a free step
                modelN.ZJ0fit = sum(((modelN.ZJ0 - ZJ0.value)./ZJ0_errorN).^2);
                modelN.fit    = modelN.ZJ0fit;
                modelN.nfit   = modelN.fit/(length(ZJ0.value));
                modelN.llh    = -0.5*(log(2*pi)*(length(modelN.ZJ0)) + ...
                    sum(log(ZJ0_errorN)) - modelN.fit);
                %uncorrelated errors, so the determinent is just the product of
                %the trace, here as a log
                detold = sum(log(ZJ0_error));
                detnew = sum(log(ZJ0_errorN));
                
            end
                
            alpha = min([ 0 (2*(detold - detnew) - (modelN.fit - model.fit)/2)]);%2 not in the previous sums
            
            if r <= alpha
                
                model = modelN;
                
            end
            
    end
                                   
end

%for change
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


%for move
%             else %delayd rejection
%                
%                 modelN2 = model;
%                 
%                 pert                  = normrnd(0, inverse_parameters.sig.depths/5);
%                 modelN2.(mode).z(move) = modelN2.(mode).z(move) + pert;
% 
%                 if modelN2.(mode).z(move) < 0 || modelN2.(mode).z(move) > inverse_parameters.depth
% 
%                     return
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
%                 alpha = min([ 1 exp(-(modelN2.fit - model.fit)/2)*(exp( - ((5*(modelN.(mode).value(move) - modelN2.(mode).value(move)).^2)/inverse_parameters.sig.(mode))))/...
%                     (exp( - (((modelN.(mode).value(move) - model.(mode).value(move)).^2)/inverse_parameters.sig.(mode))))* ...
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
% 


%errors
%         case 5 %hyperparameters
%                                                 
%             par    = rand(1);
%             modelN = model;
%             
%             if strcmp(inverse_parameters.H.prior_distribution, 'log-uniform')
%             
%                 modelN.HV_error = log(modelN.HV_error);
%                 modelN.ZJ0_error = log(modelN.ZJ0_error);
%                 
%             end
%                 
%             if par < 0.5
%                 
%                 modelN.HV_error = max([ inverse_parameters.H.sigHV_limits(1) ...
%                     min([ inverse_parameters.H.sigHV_limits(2) ...
%                     normrnd(modelN.HV_error, inverse_parameters.H.sigHV) ]) ]);
%                 
%             elseif par >= 0.5
%                 
%                 modelN.ZJ0_error = max([ inverse_parameters.H.sigZJ0_limits(1) ...
%                     min([ inverse_parameters.H.sigZJ0_limits(2) ...
%                     normrnd(modelN.ZJ0_error, inverse_parameters.H.sigZJ0) ]) ]);
%                 
%             end
%             
%             if strcmp(inverse_parameters.H.prior_distribution, 'log-uniform')
%                 
%                 modelN.HV_error  = exp(modelN.HV_error);
%                 modelN.ZJ0_error = exp(modelN.ZJ0_error);
%                 
%             end
%             
%             modelN.HVfit  = sum(((modelN.HV - HV.value)/modelN.HV_error).^2);
%             modelN.ZJ0fit = sum((( modelN.ZJ0 - ZJ0.value)/modelN.ZJ0_error).^2);
%             
%             modelN.fit = (modelN.HVfit + modelN.ZJ0fit);
%             
%             %same in both
%             alpha = min([ log(1) (log(model.HV_error/modelN.HV_error) + log(model.ZJ0_error/modelN.ZJ0_error) + ...
%                 length(modelN.HV)*log(model.HV_error/modelN.HV_error) ...
%                 + length(modelN.ZJ0)*log(model.ZJ0_error/modelN.ZJ0_error) ...
%                 - (modelN.fit - model.fit)/2 ) ]);
%             
%             r = log(rand(1));
%             if r <= alpha
%                 
%                 model = modelN;
%                 
%             end


%         case 5 %collect/seperate adjacent nodes
%             
%             while ~modelN.valid %techincally could be infinite
%             
%                 modelN       = model; %will loop, and don't want to build up pert
%                 modelN.valid = 1;
% 
%                 change1 = randi(modelN.(mode).n, 1);
%                 
%                 pert = normrnd(0, sig);
%                 
%                 modelN.(mode).value(change1)   = modelN.(mode).value(change1) + pert;
%                 
%                 if change1 == 1 %top node
%                 
%                     change2 = 1;
%                     
%                 elseif change1 == modelN.(mode).n %bottom node
%                     
%                     change2 = -1;
%                     
%                 else
%                    
%                     change2 = (-1).^(randi([1 2] ,1));
%                     
%                 end
%                 
%                 modelN.(mode).value(change1 + change2)   = modelN.(mode).value(change1 + change2) - pert;%opposite sign
%                 
%                 if modelN.(mode).value(change1) < inverse_parameters.limits.(mode)(1) || modelN.(mode).value(change1) > inverse_parameters.limits.(mode)(2)
%                     
%                     modelN.valid = 0;
%                     
%                 end
%                 
%                 if modelN.(mode).value(change1 + change2) < inverse_parameters.limits.(mode)(1) || modelN.(mode).value(change1 + change2) > inverse_parameters.limits.(mode)(2)
%                     
%                     modelN.valid = 0;
%                     
%                 end
%                 
%                 if modelN.valid
%                 
%                     [ modelN ] = evaluate_models_matdisp(modelN, HV, ZJ0, inverse_parameters);
% 
%                 end
%                     
%                 count = count + 1;
%                 
%                 if count > max_count
%                     
%                     return
%                     
%                 end
%                 
%             end
%                 
%             alpha = min([1 exp(-(modelN.fit - model.fit)/2) ]);
%             
%             r = rand(1);
%             
%             if r <= alpha
%                 
%                 model = modelN;
%                 
%             end
