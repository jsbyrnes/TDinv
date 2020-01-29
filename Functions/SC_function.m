function [ SC ] = SC_function( Z, HN, HE, window_length, onebit, nsmooth, low, high, time_keep)
%data should be iris fetch structures from just two stations. 
    
    %make the windows first
    %convert to time 
    window_length = window_length*Z.sampleRate;
    
    if mod(window_length, 2) == 0%make it odd
        
        window_length = window_length - 1;
        
    end        
    
    win_start     = 1:window_length:( length(Z.data) - window_length);
                       
    for k = 1:length(win_start)
    
        Zwin  = Z;
        HNwin = HN;
        HEwin  = HE;
        
        Zwin.data = detrend(Zwin.data(round(win_start(k)):...
            round((win_start(k) + window_length - 1))));
        HNwin.data = detrend(HNwin.data(round(win_start(k)):...
            round((win_start(k) + window_length - 1))));
        HEwin.data = detrend(HEwin.data(round(win_start(k)):...
            round((win_start(k) + window_length - 1))));
        
        Zwin.sampleCount = length(Zwin.data);
        HNwin.sampleCount = length(HNwin.data);
        HEwin.sampleCount = length(HEwin.data);
                
        if onebit
           
            Zwin.data  = sign(Zwin.data);
            HNwin.data = sign(HNwin.data);      
            HEwin.data = sign(HEwin.data);      
                        
        end
        
        Zwin  = wfFFT2(Zwin);
        HNwin = wfFFT2(HNwin);
        HEwin = wfFFT2(HEwin);
        
        Zdenom = smooth(abs(Zwin.DFT), nsmooth, 'lowess').^2;
        Ndenom = smooth(abs(HNwin.DFT), nsmooth, 'lowess').^2;
        Edenom = smooth(abs(HEwin.DFT), nsmooth, 'lowess').^2;
                
        ZZ(:, k) = real(fftshift(ifft((Zwin.DFT.*conj(Zwin.DFT))./Zdenom)));
        ZN(:, k) = real(fftshift(ifft((HNwin.DFT.*conj(Zwin.DFT))./Zdenom)));
        ZE(:, k) = real(fftshift(ifft((HEwin.DFT.*conj(Zwin.DFT))./Zdenom)));
        NN(:, k) = real(fftshift(ifft((HNwin.DFT.*conj(HNwin.DFT))./Ndenom)));
        NZ(:, k) = real(fftshift(ifft((Zwin.DFT.*conj(HNwin.DFT))./Ndenom)));
        NE(:, k) = real(fftshift(ifft((HEwin.DFT.*conj(HNwin.DFT))./Ndenom)));
        EE(:, k) = real(fftshift(ifft((HEwin.DFT.*conj(HEwin.DFT))./Edenom)));
        EZ(:, k) = real(fftshift(ifft((Zwin.DFT.*conj(HEwin.DFT))./Edenom)));
        EN(:, k) = real(fftshift(ifft((HNwin.DFT.*conj(HEwin.DFT))./Edenom)));
        
        SC(:, k, 1, 1) = butterworthfilt( ZZ(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 1, 2) = butterworthfilt( ZN(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 1, 3) = butterworthfilt( ZE(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 2, 2) = butterworthfilt( NN(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 2, 1) = butterworthfilt( NZ(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 2, 3) = butterworthfilt( NE(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 3, 3) = butterworthfilt( EE(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 3, 1) = butterworthfilt( EZ(:, k), 1/Z.sampleRate, low, high);
        SC(:, k, 3, 2) = butterworthfilt( EN(:, k), 1/Z.sampleRate, low, high);
        
    end
            
    n        = length(ZZ) - 1;
            
    ind_keep = round(n/2) - time_keep*Z.sampleRate:(round(n/2) + time_keep*Z.sampleRate);
        
    SC = SC(ind_keep, :, :, :);
    
end

