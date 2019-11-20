function [ answer, incident_phase ] = makerf(tmpP, tmpSV, nevents, Model, zerot, tlength, rp, dscr, do_pretaper)
%makerf Function to run the source convolutions, then do the deconvolution
%and align/filter/cut

    dt = Model.dt;
    
    slen = length(tmpP);
    
    %flip the time series if it is an SV arrival
    
    if strcmpi(Model.phase, 'SV')
        
        tmpP = flipud(tmpP);
        tmpSV = flipud(tmpSV);
       
        %There's a problem here. Later, the code makes a pick for the
        %receiever fucntion, but cause an array bounds error because of
        %this time flip. This solved by a circshift here that is just a
        %largeish number, there green's functions should be mostly zeros
        %and rather insensitive to how this is done.
        
%         tmpP = circshift(tmpP, round(length(tmpP)/2));
%         tmpSV = circshift(tmpSV, round(length(tmpSV)/2));
        
    end
    
    %now zero pad
    
    if Model.use_source_estimate
        
        len = length(dscr) - length(tmpP);
        
    else
        
        len = length(tmpP);
        
    end
        
    tmpP = [tmpP; zeros(len, 1)];
    tmpSV = [tmpSV; zeros(len, 1)];
            
    FtmpP = fft(tmpP);
    FtmpSV = fft(tmpSV);
    
    P_syn = zeros(length(tmpP), nevents);
    SV_syn = zeros(length(tmpSV), nevents);
    
    %make the traces by frequency domain convolution
    
    for j = 1:nevents
        
        %use the source estimate for the real inversion
        %might want to use, might not want to use
        
        if Model.use_source_estimate
        
            P_syn(:, j) = real(ifft(dscr(:, j).*FtmpP));
            SV_syn(:, j) = real(ifft(dscr(:, j).*FtmpSV));
            
        else
            
            P_syn(:, j) = tmpP;
            SV_syn(:, j) = tmpSV;
        
        end
        
        [ P_syn(:, j), SV_syn(:, j) ] = ZR2PSV( Model.phase, P_syn(:, j), SV_syn(:, j) , rp, Model.vincident(1));
        
    end
        
    %do the prefilter
    [b,a] = butter(2,[.01 2].*(2*dt)); % filter coeffs

    for i = 1:nevents
 
       P_syn(:, i) = single( filtfilt( b,a,double(P_syn(:, i)) ) );
       SV_syn(:, i) = single( filtfilt( b,a,double(SV_syn(:, i)) ) );
    
    end
    
    %now you need to cut down the SV greens
    if strcmpi('SV', Model.phase)
        
        %redefine ind, since you've move the shit out of it at this point
        %ind = ind + slen/2;
        
        %pick per event
                
        for ii = 1:nevents
        
            %pick = ind + round(timed_picks(ii)/dt);%This is to simulate the pick that you made in med_bow
                                                    %Acts like the pick from
                                                    %Cull
            
            %pick = 2450;
            
            %try a new way to do the pick!!
            [~, ind] = max(abs(SV_syn(:, ii)));
            
            %pick = ind + 2*Model.period/Model.dt;
            pick = ind - 2*Model.period/Model.dt;
                                                    
            if pick < length(P_syn)

                try
                                                    
                    Pcut(:, ii) = P_syn(pick: pick + slen - 1, ii);
                    SVcut(:, ii) = SV_syn(pick : pick + slen - 1, ii);
                
                catch
                    
                    dummy = 1;
                    
                end

            
            else
                
                Pcut(:, ii) = [P_syn(pick - slen + 1 : end, ii); zeros(slen - length(P_syn(pick - slen + 1 : end, ii)), 1)];
                SVcut(:, ii) = [P_syn(pick - slen + 1 : end, ii); zeros(slen - length(P_syn(pick - slen + 1 : end, ii)), 1)];

            end
        
        end
                                
    end
    
    P_syn = detrend(P_syn, 'constant');
    SV_syn = detrend(SV_syn, 'constant');
    
    %taper the traces
    
    r = 5/(dt*.5*length(P_syn));
    tap = tukeywin(length(P_syn), r);
    
    tap = repmat(tap, 1, nevents);
    
    P_syn = P_syn.*tap;
    SV_syn = SV_syn.*tap;
    
    if Model.use_source_estimate
                
        if strcmpi('P', Model.phase)
            
            new_source = make_source_struct(P_syn, SV_syn, dt, Model);
            
        elseif strcmpi('SV', Model.phase)
            
            new_source = make_source_struct(SV_syn, P_syn, dt, Model);
            
        end
        
    end
    
    if Model.phase == 'P'
        
        switch Model.deconvolution_method
            
            case 1 % Mercier et al., 2006. LSQ and minimum phase assumption
                
                if Model.use_source_estimate
                    
                    [G, ~] = LSQdecon(P_syn, SV_syn, zeros(size(P_syn)), new_source, length(P_syn(:, 1)), dt, dt*zerot, Model.damping_factor);
                    
                else
                    
                    [G, ~] = LSQdecon(P_syn, SV_syn, zeros(size(P_syn)), [], length(P_syn(:, 1)), dt, dt*zerot, Model.damping_factor);
                    
                end
                
            case 2 % Helffrich, 2006. Extended Time Multitaper
                
                G = drive_mtm(P_syn, SV_syn, zeros(size(SV_syn)), Model, zerot, E, V);
                
        end
                

    elseif strcmpi(Model.phase, 'SV')
        
        switch Model.deconvolution_method
            
            case 1
                
                if Model.use_source_estimate
                    
                    [G, ~] = LSQdecon(SV_syn, P_syn, zeros(size(P_syn)), new_source, length(P_syn(:, 1)), dt, dt*zerot, Model.damping_factor);
                    
                else
                    
                    [G, ~] = LSQdecon(SV_syn, P_syn, zeros(size(P_syn)), [], length(P_syn(:, 1)), dt, dt*zerot, Model.damping_factor);
                    
                end
                
            case 2
                
                G = drive_mtm(SV_syn, P_syn, zeros(size(SV_syn)), Model, zerot);
                
        end
        
    end

    %cut out zero padding

    G.C1 = G.C1(1:end-len);
    G.C2 = G.C2(1:end-len);
    
    G.C1 = detrend(G.C1);
    G.C2 = detrend(G.C2);
    %G.C3 = detrend(G.C3);
    
    r = 5/(dt*.5*length(G.C1));
    tap = tukeywin(length(G.C1), r);
    
    G.C1 = G.C1.*tap;
    G.C2 = G.C2.*tap;
    
    %done twice
%     [amp, ~] = max(G.C1);
% 
%     G.C1 = G.C1/amp;
%     G.C2 = G.C2/amp;

    G.C1 = bandpassfilt(G.C1, dt, Model.high, Model.low);
    
    %put the taper in C2 - don't need it on the other one
    
    if do_pretaper
    
        ptap = zeros(length(G.C2), 1);
    
        ptap(150:150 + (5-dt)/dt) = .001*exp(.2*log(1000)*(0:dt:(5-dt)))';
    
        ptap(150 + 5/dt:end) = 1;

        G.C2 = G.C2.*ptap;
        
    end
            
    G.C2 = bandpassfilt(G.C2, dt, Model.high, Model.low);
    %G.C3 = bandpassfilt(G.C3, dt, high, low);
    
    [amp, ~] = max(abs(G.C1));

    G.C1 = G.C1/amp;
    G.C2 = G.C2/amp;
    %G.C3 = G.C3/amp;
    
    %now cut it down to size
    
    if strcmpi('SV', Model.phase)
        
        G.C2 = G.C2*-1;
        
    end
    
    answer = G.C2(1:tlength);
    incident_phase = G.C1(1:tlength);

end