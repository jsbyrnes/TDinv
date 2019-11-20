function [ ts, SOS, G] = butterworthfilt( ts, dt, low, high, SOS, G )
%BANDPASSFILT Bandpass filters the data between high and low in Hz
%Data is tapered first with a r = .1 Tukey window(hardwired, but easy to
%change)

    n = length(ts);
    tap = tukeywin(n, .1);

    if nargin == 4
    
        [z, p, k] = butter(8,[low high].*(2*dt));
        [SOS, G] = zp2sos(z, p, k);
    
    end
        
    ts = filtfilt(SOS, G,ts.*tap);
        
end

