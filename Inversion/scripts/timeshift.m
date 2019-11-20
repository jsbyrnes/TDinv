function [ ts ] = timeshift( ts, n, nshift )
%TIMESHIFT Shift time series in the positive direction by shift/dt samples
%   Decovolution will 'shift' time series in the negative direction if time
%   series are aligned. For visualtion purposes you do not want the zero
%   point to be at n = 1; shift is in seconds
%   At the moment, this is only implemented for ts to be a vector. Note
%   that you lose any data in the last shift seconds of the time series
    
    ts = [zeros(nshift, 1); ts(1:n - nshift)]';

end