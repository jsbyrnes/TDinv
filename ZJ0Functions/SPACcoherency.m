function [ f, Z1, Z2 ] = SPACcoherency( Z1, H1, H2, ...
    window_length, normalize, taper_percent)
%data should be iris fetch structures from just two stations. 
    
    %make the windows first
    %convert to time
    window_length = window_length*Z1.sampleRate;
    win_start     = 1:window_length*(1 - overlap):( length(Z1.data) - window_length);
    taper         = tukeywin(window_length, taper_percent);%hanning window at 1
        
    if integrate
       
        Z1.data = cumsum(Z1.data);
        H1.data = cumsum(H1.data);
        H2.data = cumsum(H2.data);
        
    end
    
    if normalize
        
        Z1.data = Z1.data./abs(hilbert(Z1.data));
        H1.data = H1.data./abs(hilbert(H1.data));
        H2.data = H2.data./abs(hilbert(H2.data));
        
    end
       
    %make the window
    
    for k = 1:length(win_start)
    
        Z1win = Z1;
        Z2win = Z2;
        Z1win.data = detrend(Z1win.data(round(win_start(k)):...
            round((win_start(k) + window_length - 1))));
        Z2win.data = detrend(Z2win.data(round(win_start(k)):...
            round((win_start(k) + window_length - 1))));
        Z1win.data = taper.*Z1win.data;
        Z2win.data = taper.*Z2win.data;        
        Z1win.sampleCount = length(Z1win.data);
        Z2win.sampleCount = length(Z2win.data);
                
        Z1win = wfFFT2(Z1win);
        Z2win = wfFFT2(Z2win);
        
        CZall(:, k) = (Z1win.DFT.*conj(Z2win.DFT))./sqrt(Z1win.DFT.*conj(Z1win.DFT)...
            .*Z2win.DFT.*conj(Z2win.DFT));
            
    end
        
    %CZ = nanmean(atanh(CZ), 2);
    %CZ       = tanh(nanmean(atanh(real(CZall)), 2)) + 1i*tanh(nanmean(atanh(imag(CZall)), 2));
    %CZ_error = std(bootstrp(1e4, @(x) tanh(median(atanh(real(x)))), CZall));
   
    CZ       = tanh(nanmedian(real(atanh(real(CZall))),2)) + 1i*nanmedian(imag(CZall),2);        
    %get error by bootstrapping
    
    for k = 1:length(CZ)
    
        CZ_error(k) = std(bootstrp(1e2, @(x) tanh(nanmedian(real(atanh(real(x))))), CZall(k, :)));
    
    end
        
    if ~isempty(H1)
    
        R1.data = detrend(R1.data);
        R2.data = detrend(R2.data);
        T1.data = detrend(T1.data);
        T2.data = detrend(T2.data);

        if normalize

            R1.data = R1.data./abs(hilbert(R1.data));
            R2.data = R2.data./abs(hilbert(R2.data));
            T1.data = T1.data./abs(hilbert(T1.data));
            T2.data = T2.data./abs(hilbert(T2.data));

        end
                       
        R1 = wfFFT2(R1);
        R2 = wfFFT2(R2);
        CR = (R1.DFT.*conj(R2.DFT))./sqrt(R1.DFT.*conj(R1.DFT)...
            .*R2.DFT.*conj(R2.DFT));
        
        T1 = wfFFT2(T1);
        T2 = wfFFT2(T2);
        CT = (T1.DFT.*conj(T2.DFT))./sqrt(T1.DFT.*conj(T1.DFT)...
            .*T2.DFT.*conj(T2.DFT));
        
    else
        
        CR = [];
        CT = [];
        
    end
       
    f = Z1win.frequencies;
    
end

