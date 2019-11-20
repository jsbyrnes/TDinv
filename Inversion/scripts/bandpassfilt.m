function [ ts ] = bandpassfilt( ts, dt, high, low )
%BANDPASSFILT bandpassfilt( ts, dt, high, low )
%Bandpass filters the data between high and low in Hz
%Data is tapered first with a r = .1 Tukey window(hardwired, but easy to
%change). Assumes data is in the columns!!

%Joseph Byrnes jbyrnes@uoregon.edu July 2012

    [n q] = size(ts);
    
    r = 5/(dt*.5*(n));
    
    tap = tukeywin(n, r);
    
    [b,a] = butter(2,[low high].*(2*dt));

    for i = 1:q    
        
        ts(:, i) = single(filtfilt(b,a,double(ts(:, i).*tap).*tap));
        
    end
    
end

