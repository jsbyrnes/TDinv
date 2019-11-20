function [ Pcon, SVcon , source] = convolvesource( Pgrns, SVgrns, source, phase)
%CONVOLVESOURCE Convolve the green's function with a given source or with
%the simple source that has been being used
%Green's functions should be columns

    [m, q] = size(Pgrns);

    n = length(source);
    
    if n > m
        
        extra = n - m;
        
        Pgrns = [Pgrns; zeros(extra,q)];
        SVgrns = [SVgrns; zeros(extra,q)];
        
    elseif m > n
        
        extra = m - n;
        
        Pgrns(m - extra, :) = [];
        SVgrns(m - extra, :) = [];
        
    end
    
    %I want to zero pad but I can't change the length so circshift
    %These green's functions are all the same
    
    if strcmpi('P', phase)
        
        [~, zeropoint] = max(Pgrns(:, 1));
        
    elseif strcmpi('SV', phase)
        
        [~, zeropoint] = max(SVgrns(:, 1));

    end
    
    Pgrns = circshift(Pgrns, [(-zeropoint + 1) 0]);
    SVgrns = circshift(SVgrns, [(-zeropoint + 1) 0]);

    F_P = fft(Pgrns, [], 1);
    F_SV = fft(SVgrns, [], 1);
    
    if isempty(source)
        
        load('source45966.mat');
       
        source = source45966;
        
    end
        
    F_source = fft(source);
    
    for i = 1:q
    
        Pcon(:, i) = F_P(:,i) .* F_source;
        SVcon(:, i) = F_SV(:,i) .* F_source;
        
    end
    
    Pcon = real(ifft(Pcon, [], 1));
    SVcon = real(ifft(SVcon, [], 1));

end