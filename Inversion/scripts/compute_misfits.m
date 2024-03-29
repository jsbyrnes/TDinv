function [ RF, Greens ] = compute_misfits( RF, Greens, high, low, Model, plotflag, saveplots, sourcename)
%COMPUTE_MISFIT Fills in the rms vectors by repeated
%convolution/deconvolution of the filtered synthetics.

    rn = length(RF);
    gn = length(Greens);

    load(sourcename);
    
    %I need this to help out the parfor loop
    %can't classify the Source variable - needs to know if it's a function
    %or a variable
    %Sources = Source;
    
    %to use parfor, you have to avoid using the Model structure
    %not sure where it was erroring, so I just removed (almost) all of it
    %from the for loop
    dt = Model.dt;
    phase = Model.phase;
    rp = Model.rp;
    len = Model.len;
    ini_range = Model.range;
    
    source_len = length(Source(1).td_source_est);
    
    %do the same thing for some green's structure vectors
    greens_rp = [Greens(:).rp];

    %I don't know if this is the best way to do this, feel free to change
    %it
    %temp array to help index the source struct
    
    for i = 1:length(Source)

        eid(i) = Source(i).eid; %#ok<AGROW>

    end

    for ri = 1:1
        
        if RF(ri).exists == 0
            
            continue
            
        end

        fprintf('Running Receiver Function %i\n', ri);
        
        %same as with model
        events = RF(ri).events;
        rf_zerot = RF(ri).zerot;
        rf_trace = RF(ri).trace;
        
        rf_rms = zeros(1, gn);
        t0 = zeros(1, gn);
        
        %first, find the index of the correct ray parameter in the model
        %structure
        
        [~, rf_rp_index, ~ ] = findnearest(rp, RF(ri).rp, []);
        
        tmp_obs = cell(gn, 1);
        tmp_syn = cell(gn, 1);
        
        nevents = length(events);
        
        ip = zeros(1, nevents);
                            
        for j = 1:nevents
            
            ip(j) = find(eid == events(j));
            
        end
        
        %only need one set of sources for each rf
        source_cut = Source(ip);
        
        h = waitbar(0, 'Running current rf');
        
        parfor gi = 1:gn
                                                
            %very first thing to do is check if the rp is right by seeing
            %if it is the clostest one in the structure
            %don't need to consider any others
            
            [~, gr_rp_index, ~ ] = findnearest(rp, greens_rp(gi), []);
            
            if rf_rp_index ~= gr_rp_index
                
                %This means we should not even consider this one - no point
                
                tmp_obs{gi} = 0;
                tmp_syn{gi} = 0;

                Greens(gi).rms(ri) = NaN;
                rf_rms(gi) = NaN;
                
                continue
                
            end
            
            %make the greens function now, save as local greens(don't save the traces after the loop, takes lots of ram)
            
            [z, vp, vpvs] = compile_velocity_model(Greens(gi));
            
            %make rho with the brocher curve, for now just the crappy one
            rho = .33 + .77*vp;
            
            [syn_P, syn_SV, ~] = anirec(phase, dt, greens_rp(gi), 0, vp, vp./vpvs, z, rho, zeros(1, length(vp)), ...
                zeros(1, length(vp)), zeros(1, length(vp)), zeros(1, length(vp)), zeros(1, length(vp)));
            
            P_syn = zeros(source_len, nevents);
            SV_syn = zeros(source_len, nevents);

            for j = 1:nevents

                [P_syn(:, j), SV_syn(:, j), ~] = convolvesource(syn_P, syn_SV, source_cut(j).td_source_est, phase);

            end
                 
            if phase == 'P'
                
                [G, ~] = LSQdecon(P_syn, SV_syn, zeros(size(P_syn)), source_cut, length(P_syn(:, 1)), dt, dt*400);
                
            elseif strcmpi(phase, 'SV')
                
                [G, ~] = LSQdecon(SV_syn, P_syn, zeros(size(P_syn)), source_cut, length(P_syn(:, 1)), dt, dt*400);
                
            end
                
            %[G, ~] = LSQdecon(P, SV, zeros(size(P)), RF(ri).source, length(P(:, 1)), Model.dt, 0);

            %[G, ~] = LSQdecon(P, SV, zeros(size(P)), [], Model.len, Model.dt, Model.dt*RF(ri).zerot);

            %for some reason, LSQDecon shifts the G.C2/3 system too much
            %I wish it didn't shift it at all ever but oh well
            
            %Idk if this is necessary but why the hell not
            %Probably should is using RFs from medb_bow since this step is
            %applied to the observed ones
            G.C1 = detrend(G.C1);
            G.C2 = detrend(G.C2);
            %G.C3 = detrend(G.C3);

            G.C1 = bandpassfilt(G.C1, dt, high, low);
            G.C2 = bandpassfilt(G.C2, dt, high, low);
            %G.C3 = bandpassfilt(G.C3, dt, high, low);

            %need to rezero and normalize around 400 - the hardwired point
            %that the pulse should have arrived at
            %give it a tolerence of ~2 seconds, doesn't need to be exact
            [amp, syn_zero] = max(G.C1(400 - 1/dt:400 + 1/dt));
            
            %should come to exactly 400 but whatever I already wrote it and
            %I am a lazy bastard
            syn_zero = syn_zero + (400 - 1/dt - 1);
            
            Greens(gi).zerot = syn_zero;

            G.C1 = G.C1/amp;
            G.C2 = G.C2/amp;
            %G.C3 = G.C3/amp;
            
            %clip synthetic to be the same length as a real one
            %might need to zero pad the left hand side if too short
            %all that's important that the zeros are the same index
                        
            if isempty(ini_range)
                %range is in seconds. If you make the range empty in
                %define_project, the code will make it do the misfit over
                %the rf shift to the end of either the receiver functions
                %or the synthetic, which ever is shorter. Assumes the S
                %receiver function is in negative time and the S green's
                %function is in real time

                %Use Model.dt to convert these to seconds. This makes it
                %easier
                if phase == 'P'

                    a1 = dt*max([ -abs(rf_zerot) -abs(Greens(gi).zerot)]);
                    a2 = dt*min([length(rf_trace(rf_zerot:end)) length(Greens(gi).P(Greens(gi).zerot:end))]); %#ok<*PFBNS>

                    range = [ a1 a2 ];

                elseif strcmpi(phase, 'SV')

                    a1 = -dt*min([ abs(rf_zerot) ]);
                    a2 = dt*min([length(rf_trace(rf_zerot:end)) length(Greens(gi).P(1:Greens(gi).zerot))]);

                    range = [ a1 a2 ];

                end

            elseif (abs(ini_range(1)) > dt*rf_zerot) || (abs(ini_range(1)) > dt*Greens(gi).zerot)

                    range(1) = dt*max([ -abs(rf_zerot) -abs(Greens(gi).zerot)]);
                    range(2) = ini_range(2);

            elseif (abs(ini_range(2)) + dt*rf_zerot >= len) || (abs(ini_range(2)) + dt*Greens(gi).zerot >= len)

%                     Model.range(2) =

            end

            SynRF = G.C2;

            ObsRF = rf_trace;

            %get the section of the RF and syn RF that we need

            bound1 = round(abs(range(1)/dt));

            bound2 = round(abs(range(2)/dt));

            %Important, this comment is a repeat but it's important. Code
            %assumes that the ObsRF is negative time(flipped), but the
            %synthetic is real time(not flipped) for SRF;

            %first, get the right section of the observed RF

            obs_sec = ObsRF(rf_zerot - bound1 + 1 : rf_zerot + bound2 - 1);

            %now the synthetic section, which is different for P and S

            if phase ==  'P'

                syn_sec = SynRF(Greens(gi).zerot - bound1 + 1 : Greens(gi).zerot + bound2 - 1);

            elseif strcmpi(phase, 'SV')

                %flip it to match the ObsRF
                syn_sec = wrev(SynRF(Greens(gi).zerot - bound2 + 1: Greens(gi).zerot + bound1 - 1));
                %syn_sec = SynRF(Greens(gi).zerot - bound2 + 1: Greens(gi).zerot + bound1);

            end

            %obs_sec and syn_sec should be align perfectly in time.

            diff = syn_sec-obs_sec;
            
            %apply an envelope
            
            envelope = zeros(length(diff), 1);
            
            if strcmpi(phase, 'P')
            
                envelope(bound1 - .5/dt:bound1 + 3/dt) = tukeywin(length(bound1 - .5/dt:bound1 + 3/dt), .1);
                
            elseif strcmpi(phase, 'SV')
                
                envelope(bound2 - .5/dt:end) = 1;
                
            end
            
            diff = diff.*envelope;

            rms = sqrt(mean((diff.^2)));

            %Save the RMS misfits in two different place, both indexed
            %in a way that's pretty easy to access
            rf_rms(gi) = rms;
            t0(gi) = bound1;
            Greens(gi).rms(ri) = rms;

            tmp_obs{gi} = obs_sec;
            tmp_syn{gi} = syn_sec;
            
            waitbar(gi/gn, h);

        end
        
        delete(h)
        
        RF(ri).rms = rf_rms;
        
        %make these a value way too large to even happen
        RF(ri).rms(isnan(RF(ri).rms)) = 1;
        
        [~, ind] = min(RF(ri).rms);

        RF(ri).bestfit_obs = cell2mat(tmp_obs(ind));
        RF(ri).bestfit_syn = cell2mat(tmp_syn(ind));
        
        if plotflag == 'y'

            
            %%%%%%%
            %figure 23 compares the best fitting synthetic with the
            %observed receiver function
            %%%%%%%
            figure(23)
            clf
            
            len = length(RF(ri).bestfit_obs);
            t = -((t0(ind) - 1)*dt):dt:(len-t0(ind))*dt;
            
            envelope = zeros(1,len);
            
            envelope(t0(ind) - .5/dt:len) = tukeywin(length(t0(ind) - .5/dt:len), .1);
            
            %normalize the window being used here to the max of the rfs so
            %that you can actually read the graph that's being spat out at
            %you
            
            norm = max([RF(ri).bestfit_obs; RF(ri).bestfit_syn]);
            
            envelope = envelope*norm;
            
            plot(t, RF(ri).bestfit_obs, 'k');
            hold on
            plot(t, RF(ri).bestfit_syn, 'r');
            plot(t, envelope, 'b--');
            xlabel('Time(seconds)')
            title(['RMS = ' num2str(RF(ri).rms(ind)) ' and Greens is ' num2str(ind)]);
            
            %%%%%%%%%%
            %figure 24 gives the velocity model associated with the best
            %fitting synthetic
            %%%%%%%%%%
            figure(24)
            clf
            
            [z, vp, vpvs] = compile_velocity_model(Greens(ind));
            
            z(length(z)) = 100;
            
            z = [ 0 z 100];
            
            vp = [ vp(1) vp  vp(length(vp))];
            
            vpvs = [ vpvs(1) vpvs vpvs(length(vpvs))];
            
            vs = vp./vpvs;
            
            stairs(vp, z, 'b-', 'LineWidth', 2), hold on, stairs(vs, z, 'r-', 'LineWidth', 2)
            
            set(gca, 'YDir', 'reverse'), grid on
            
            xlim([(min(vs) - 1) (max(vp) + 1)]);
            
            ylim([ 0 z(length(z)-1)/2 ]);
            
            title('Choosen Velocity Model'), ylabel('Depth(Kilometers)'), xlabel('Velocity(Km/s)')
            
            %%%%%%%%%%
            %figure 25 plots out the rms space with a histogram
            %pretty basic, will need something fancier in the future
            %%%%%%%%%%
            
            figure(25)
                        
            %first get the values that aren't 1, these are the points that
            %are skipped by this rf because the ray parameter was wrong
            
            rmsvector = RF(ri).rms(RF(ri).rms ~= 1);
            numberofgreens = length(rmsvector);
            
            numbins = numberofgreens/100;
            
            if numbins < 10
                
                numbins = 10;
                
            end
            
            hist(rmsvector, numbins)
            
            xlabel('RMS values'), title(['Histogram of RMS values out of ' num2str(numberofgreens) ' tested'])
                        
            drawnow

        end
       
        if plotflag =='y' && saveplots == 'y'
           
            print(23, ['Plots/rf#' num2str(ri) '.bestfitresult.tiff'], '-dtiff');
            print(24, ['Plots/rf#' num2str(ri) '.bestfitmodel.tiff'], '-dtiff');
            print(25, ['Plots/rf#' num2str(ri) '.rmshistogram.tiff'], '-dtiff');
            
        end

    end

end
