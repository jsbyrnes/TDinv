function [ C, C_error ] = SPAC_butterworth( data1, data2, sample_rate, segment_length, ...
    filter_width, filters, normalize, make_complex, integrate, onebit, bts)
%data should be iris fetch structures from just two stations. 
        
    for i = 1:(length(filters))
                                             
        %50% overlap
        n_inc     = floor(length(data1)/(segment_length*sample_rate/filters(i))) - 1;
        n_samples = floor(segment_length*sample_rate/filters(i));
        
        data1_seg = data1( 1:n_samples);
        [~, SOS, G] = butterworthfilt(data1_seg, 1/sample_rate, filters(i) - filter_width(i)/2, filters(i) + filter_width(i)/2);
        
        for j = 1:n_inc

            %%%%%%%%%%%%%
            %pair 1 and 2
            seg1 = data1((j-1)*n_samples+1:j*n_samples);
            seg2 = data2((j-1)*n_samples+1:j*n_samples);
                        
            seg1 = detrend(seg1);
            seg2 = detrend(seg2);
            
            if integrate
               
                seg1 = cumsum(seg1);
                seg2 = cumsum(seg2);
                
            end
            
            taper = tukeywin(n_samples, 1);
                     
            seg1 = butterworthfilt(seg1.*taper,1/sample_rate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);
            seg2 = butterworthfilt(seg2.*taper,1/sample_rate, filters(i) - filter_width/2, filters(i) + filter_width/2, SOS, G);
            
            if normalize
               
                seg1 = seg1/rms(seg1);%sign(Z1_seg);%Z1_seg./abs(hilbert(Z1_seg));
                seg2 = seg2/rms(seg2);%sign(Z2_seg);%Z2_seg./abs(hilbert(Z2_seg));
                
            end
            
            if make_complex
               
                seg1 = hilbert(seg1);
                seg2 = hilbert(seg2);
                
            end
                        
            if onebit
                
                seg1 = sign(seg1);
                seg2 = sign(seg2);
                
            end
            
            cc = corrcoef(seg1,seg2);
            Call(j) = cc(1,2);
            
        end
                
        C(i)       = tanh(median(atanh(real(Call)))) + 1i*median(imag(Call));
        
        %get error by bootstrapping
        C_error(i) = std(bootstrp(bts, @(x) tanh(median(atanh(real(x)))), Call));
                 
        Call = [];
        
    end
        
end

        %pdf       = ksdensity(CZall, xi, 'BandWidth', 0.01);
        %[~, ind]  = max(pdf);
        %CZ(i)     = xi(ind); 


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

    %     CZ = mean(CZall, 2);
%     CZ_error = std(CZall, 0, 2)/sqrt(n_inc);
