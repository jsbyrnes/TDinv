function [ RF ] = populate_HV_struct( path, Model)
%POPULATE_RF_STRUCT Loads the Receiver Functions
%Path is to the file that contains the recierver functions. Right now, this
%is only implemented for the seis.mat file that comes out of med_bow. 
%Sum style is a flag that controls smoothing. 1 makes one structure for
%each RF. 2 sums together all the RFs at a given station blindly( e.g. for
%looking at shallow crustal features where the ray parameter doesn't really
%matter and the move out isn't large enough to sample much hetrogenetity).

    load(path);
    
    dt = Model.dt;
    high = Model.high;
    low = Model.low;
    pretaper = Model.do_pretaper;
    sumstyle = Model.sumstyle;
        
    [seislen, ininum] = size(seis);
    
    if Model.len == 0
        
        tracelength = seislen;
        
    else
        
        tracelength = Model.len/dt;
        
    end
    
    if pretaper
        
        %make the taper
        %looks exactly like the envelope that I've define in make_rfs
        
            ptap = zeros(seislen, ininum);
        
            ptap(150:150 + (5-dt)/dt, :) = repmat(.001*exp(.2*log(1000)*(0:dt:(5-dt)))', 1, ininum);
            
            ptap(150 + 5/dt:end, :) = 1;
            
            seis = seis.*ptap;
            
            seis = bandpassfilt(seis, dt, high, low);
            seis1 = bandpassfilt(seis1,dt, high, low);
            
            %renormalize, not sure if I have too but it's a good idea
            %note that the seis1 matrix IS NOT TAPERED
            
            amp = repmat(max(seis1, [], 1), seislen, 1);
            
            seis = seis./amp;
            seis1 = seis1./amp;
            
    end
            
    switch sumstyle
        
        case 0
            
            traces = seis(1:tracelength, :);
            traces1 = seis1(1:tracelength, :);
            
            %this gives the index vector so that it just maps one to one
            indexes = 1:ininum;
            
        case 1
            
            msta = max(istn);
            
            for i = 1:msta
                
                traces(1:tracelength, i) = sum(seis(1:tracelength, istn==i), 2);
                
                traces1(1:tracelength, i) = sum(seis1(1:tracelength, istn==i), 2);
                
                a = find(istn == i);
                
                if ~isempty(a)
                
                    indexes(i) = a(1); %this finds the meta data at the correct points
                
                elseif isempty(a)
                    
                    indexes(i) = -1; %flag that there are no measuremtns at this station
                    
                end
                
                traces(:, i) = traces(:, i)/length(a);
                traces1(:, i) = traces1(:, i)/length(a);
                
            end
            
    end
            
    [~ , numberofrfs] = size(traces);
        
    RF(1:numberofrfs, 1) = struct('trace', zeros(tracelength, 1), 'zerot', 0, 'rms', [], 'source', [], ...
        'rp', 0, 'stalat', 0, 'stalon', 0, 'olat', 0, 'olon', 0, 'sta_id', 0, 'baz', 0, 'delta', 0, 'dt', ...
        dt, 'events', [], 'dscr', []);
    
    for i = 1:numberofrfs
        
        %rms defined later
        
        ind = indexes(i);
        
        if ind ~= -1
            
            RF(i).trace = traces(:, i);
            RF(i).trace1 = traces1(:, i);
            %position, not time
            [~, RF(i).zerot] = max(traces1(:, i));
            %RF(i).zerot = 150;
            RF(i).rp = slow(ind);
            RF(i).stalat = slat(ind);
            RF(i).stalon = slon(ind);
            RF(i).olon = olon(ind);
            RF(i).olat = olat(ind);
            RF(i).sta_sid = istn(ind);
            RF(i).baz = baz(ind);
            RF(i).delta = delta(ind);
            RF(i).events = events{ind}; %this isn't actually very exact if there's a stack of events (sumstyle 2)...
            
            if Model.use_source_estimate
            
                RF(i).dscr = deconsources(ind);
            
            end
            
            RF(i).exists = 1;
            
        elseif ind == -1
            
            RF(i).exists = 0;
            
        end
        
    end
    
end

