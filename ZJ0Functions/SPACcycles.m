function [ filters, CZ, CR, CT, CZ_error, CR_error, CT_error, arclen, sta_azi ] = ...
    SPACcycles( Z1, h11, h21, Z2, h12, h22, segment_length, filter_width, ...
    freq_range, normalize, make_complex, integrate, onebit)
%data should be iris fetch structures from just two stations. 
    
    [arclen, sta_azi] = distance(Z1.latitude, Z1.longitude, Z2.latitude, Z2.longitude);
    arclen = arclen * 111.12 * 1000;
    
    %radial for one pair and not the other
    
    if ~isempty(h11)
    
        R1 = cosd(sta_azi)*h11.data + sind(sta_azi)*h21.data;
        R2 = cosd(sta_azi)*h12.data + sind(sta_azi)*h22.data;
        T1 = -1*sind(sta_azi)*h11.data + cosd(sta_azi)*h21.data;
        T2 = -1*sind(sta_azi)*h12.data + cosd(sta_azi)*h22.data;
        
    end
    
    filters = freq_range(1):filter_width:freq_range(2);%try using lots of overlap
    
    for i = 1:(length(filters))
            
        %50% overlap
        n_inc     = floor(length(Z1.data)/(segment_length*Z1.sampleRate/filters(i))) - 1;
        n_samples = floor(segment_length*Z1.sampleRate/filters(i));
           
        %to make filter
        Z1_seg = Z1.data( 1:n_samples);
        [~, SOS, G] = butterworthfilt(Z1_seg, 1/Z1.sampleRate, filters(i) - filter_width/2, filters(i) + filter_width/2);
        
        for j = 1:n_inc

            %%%%%%%%%%%%%
            %pair 1 and 2
            Z1_seg = detrend(Z1.data( (j-1)*n_samples+1:j*n_samples));
            Z2_seg = detrend(Z2.data( (j-1)*n_samples+1:j*n_samples));
            
            taper = tukeywin(n_samples, 1);
                        
            if integrate
               
                Z1_seg = cumsum(Z1_seg);
                Z2_seg = cumsum(Z2_seg);
                
            end
            
            Z1_seg = butterworthfilt(Z1_seg.*taper,1/Z1.sampleRate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);
            Z2_seg = butterworthfilt(Z2_seg.*taper,1/Z1.sampleRate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);
            
            %extra safe!
            Z1_seg = detrend(Z1_seg);
            Z2_seg = detrend(Z2_seg);
            
            if ~isempty(h11)

                R1_seg = R1.data( (j-1)*n_samples+1:j*n_samples);
                R2_seg = R2.data( (j-1)*n_samples+1:j*n_samples);
                T1_seg = T1.data( (j-1)*n_samples+1:j*n_samples);
                T2_seg = T2.data( (j-1)*n_samples+1:j*n_samples);
                
                fR1 = butterworthfilt(R1_seg,1/Z1.sampleRate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);
                fR2 = butterworthfilt(R2_seg,1/Z1.sampleRate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);
                fT1 = butterworthfilt(T1_seg,1/Z1.sampleRate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);
                fT2 = butterworthfilt(T2_seg,1/Z1.sampleRate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);

            end
            
            if normalize
               
                Z1_seg = Z1_seg/rms(Z1_seg);%sign(Z1_seg);%Z1_seg./abs(hilbert(Z1_seg));
                Z2_seg = Z2_seg/rms(Z2_seg);%sign(Z2_seg);%Z2_seg./abs(hilbert(Z2_seg));
                
            end
            
            if make_complex
               
                Z1_seg = hilbert(Z1_seg);
                Z2_seg = hilbert(Z2_seg);
                
            end
                        
            if onebit
                
                Z1_seg = sign(Z1_seg);
                Z2_seg = sign(Z2_seg);
                
            end
            
            cc = corrcoef(Z1_seg,Z2_seg);
            
            [ccc,lags] = xcorr(Z1_seg,Z2_seg);
            
            [~, ind] = max(ccc);
            
            dt(j) = lags(ind)/100;
                        
            CZall(j) = cc(1,2);

            if ~isempty(h11)
            
                R1_seg = fR1( (j-1)*n_samples+1:j*n_samples);
                R2_seg = fR2( (j-1)*n_samples+1:j*n_samples);
                cc = corrcoef(R1_seg,R2_seg);

                CRall(j) = cc(1,2);

                T1_seg = fT1( (j-1)*n_samples+1:j*n_samples);
                T2_seg = fT2( (j-1)*n_samples+1:j*n_samples);
                cc = corrcoef(T1_seg,T2_seg);

                CTall(j) = cc(1,2);
                        
            end
            
        end
        
        %pdf       = ksdensity(CZall, xi, 'BandWidth', 0.01);
        %[~, ind]  = max(pdf);
        %CZ(i)     = xi(ind); 
        
        CZ(i)       = tanh(median(atanh(real(CZall)))) + 1i*median(imag(CZall));
        
        %get error by bootstrapping
        CZ_error(i) = std(bootstrp(1e4, @(x) tanh(median(atanh(real(x)))), CZall));
        
        %CZ_error(i) = std(atanh(CZall))/sqrt(length(CZall));
        %CZ(i) = (max(pdf(round(length(pdf)/2):end) - max(pdf(1:round(length(pdf)/2)))))/(max(pdf(round(length(pdf)/2):end) + max(pdf(1:round(length(pdf)/2)))));
        %CZ(i) = mean(CZall);
        
%         for k = 1:nbs
%             
%             CZtmp = CZall(randi(n_inc, 1, n_inc));
%             pdf = ksdensity(CZall, xi, 'BandWidth', 0.01);
%             [~, ind] = max(pdf);
%             CZset(k) = xi(ind);
%             
%         %end
%         
%         %CZ_error(i) = std(CZset);
         
        CZall = [];
        CRall = [];
        CTall = [];
        
    end
    
    %     CZ = mean(CZall, 2);
%     CZ_error = std(CZall, 0, 2)/sqrt(n_inc);
    
    if ~isempty(h11)
    
        CR = mean(CRall, 2);
        CT = mean(CTall, 2);
        CR_error = std(CR, 0, 2)/sqrt(n_inc);
        CT_error = std(CT, 0, 2)/sqrt(n_inc);
    
    end
    
    if isempty(h11)
       
        CR = [];
        CT = [];
        CR_error = [];
        CT_error = [];
        
    end
    
end

