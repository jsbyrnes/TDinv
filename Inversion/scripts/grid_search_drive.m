 function [ syn_RF, model_vectors, s_rms] = grid_search_drive( RF, Model)
%grid_search_drive Performs the grid search

    %two steps to a receiver function
    %first, you generate all the Green's functions that you will need
    
    %then do the comparisons
            
    %compiles the vectors that get input into the grid search
    model_vectors = make_models(Model);
    
    total = length(model_vectors);
    
    rfn = size(RF);
                
    for ri = 1:rfn
        
        fprintf('Running RF Number %i\n', ri);
        
        %cut the rf in time
        
        answer = RF(ri).trace;
        rf_zerot = RF(ri).zerot;
        
        syn_RF = zeros(total, length(answer));
        
        dscr = []; %source from Med_bow
                   %needs to be cleared to here for the code to work
                   
        if Model.use_source_estimate
        
            dscr = RF(ri).dscr{1};
        
        end
        
%         samp = max(dscr);
%         samp = repmat(samp, length(dscr), 1);
%         
%         dscr = dscr./samp; %no reason to normalize since you can just
%         normalize after the deconvolution
        
        %make the envelope
        %right now it's hardwired, which sorta sucks but what.ev.a. for now
        
        envelope = zeros(length(RF(ri).trace), 1);
        
        if strcmpi(Model.phase, 'P')
            
            envelope(rf_zerot - .5/Model.dt:bound1 + 3/dt) = tukeywin(length(rf_zerot - .5/dt:rf_zerot + 3/dt), .1);
            
        elseif strcmpi(Model.phase, 'SV')
            
            %             envelope(rf_zerot:rf_zerot + (5-Model.dt)/Model.dt) = .001*exp(.2*log(1000)*(0:Model.dt:(5-Model.dt)));
            %
            %             envelope(rf_zerot + 5/Model.dt:end) = 1;
            
            envelope(:) = 1;
            
        end
        
        %find the sources you need to do the inversion
                
%         ip = [];
%         
%         for j = 1:nevents
%             
%             ip(j) = find(RF(ri).events(j) == [Source(:).eid]);
%             
%         end
        
        number_events = length(RF(ri).events);
        
        %generate the green's functions
        
        for j = 1:total
            
            tmp_length = length(model_vectors(j).z);
            
            zero_vec = zeros(tmp_length, 1);
            
            if strcmpi(Model.phase, 'P');
            
                rho = nafedrake_rho(model_vectors(j).vincident);
            
                [tmpP, tmpSV] = anirec(Model.phase, Model.dt, RF(ri).rp, ... 
                    RF(ri).baz, model_vectors(j).vincident,...
                    model_vectors(j).vincident./model_vectors(j).vpvs, ...
                    model_vectors(j).z, rho, zero_vec, zero_vec, zero_vec, zero_vec, zero_vec, 0);
               
                [syn_RF(j, :), syn_C1] = makerf(tmpP, tmpSV, number_events, Model, RF(ri).zerot, length(answer), RF(ri).rp, dscr, Model.do_pretaper);
            
            elseif strcmpi(Model.phase, 'SV')
               
                rho = nafedrake_rho(model_vectors(j).vincident.*model_vectors(j).vpvs);

                [tmpP, tmpSV] = anirec(Model.phase, Model.dt, RF(ri).rp, ... 
                    RF(ri).baz, model_vectors(j).vincident.*model_vectors(j).vpvs,...
                    model_vectors(j).vincident, ...
                    model_vectors(j).z, rho, zero_vec, zero_vec, zero_vec, zero_vec, zero_vec, 0);
                
                [syn_RF(j, :), syn_C1(j, :)] = makerf(tmpP, tmpSV, number_events, Model, RF(ri).zerot, length(answer), RF(ri).rp, dscr, Model.do_pretaper);
                
                dummy = 1;%just a good spot to put a breakpoint
                
            end
            
        end
                
        %grade the results
        s_rms = evaluation(syn_RF', answer, envelope)';
        
        [~, best] = min(s_rms);
        
        %plot out all the results and save everything
        
        shift_index = RF(ri).zerot;
        
        t = -shift_index*Model.dt:Model.dt:(length(RF(ri).trace)-shift_index - 1)*Model.dt;
        
        close all
        
        figure(4), clf, hold on,plot(t, syn_RF', 'Color', [.8 .8 .8]), plot(t, RF(ri).trace, 'k'), plot(t, syn_RF(best, :), 'r', 'LineWidth', 1), grid on
        plot(t, envelope*max(RF(ri).trace+.1), 'b--'); ylim([-max(RF(ri).trace + .2) max(RF(ri).trace + .2)])
        xlabel('Seconds'), title('Black is RF, red is average of the gray solutions')
        print(4, ['rfsolutions' num2str(ri) '.tiff'] , '-dtiff');
        
        save(['rfsolution' num2str(ri) '.mat'], 'model_vectors', 's_rms', 'syn_RF', 'syn_C1');
            
    end
end

