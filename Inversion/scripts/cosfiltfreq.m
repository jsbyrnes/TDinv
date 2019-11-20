function [ fs ] = cosfiltfreq( window, fs, fc, dt)
%COSFILT cosine squared filter with a cut off frequency of fc

    F = fftshift(fs);
    
    %f = (sps/2)*linspace((-1),(1),window);
    
    f = make_freq(window, dt);
    
    ff = fftshift(f);
    
    figure(5)
    
    semilogy(ff, abs(F));
    xlabel('Frequency');
    title('Unfiltered Spectra');
    
    cosfilt = zeros(window, 1);
    
    index = find(abs(ff) <= fc);
    
    cosfilt(index) = (cos(pi*ff(index)/(2*fc)).^2)/(fc*dt)';
    
    F = F.*(cosfilt);
    
    figure(6)
    semilogy(ff, abs(F));
    xlabel('Frequency');
    title('Filtered Spectra with a cutoff of fc');
    
    fs = ifftshift(F);

end

