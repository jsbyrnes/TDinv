function [ C, C_error ] = SPAC_gaussian( data1, data2, sample_rate, segment_length, ...
    filter_width, filters, normalize, make_complex, integrate, onebit, bts)
%data should be iris fetch structures from just two stations. 
        
    for i = 1:(length(filters))
                                             
        %50% overlap
        n_inc     = floor(length(data1)/(segment_length*sample_rate/filters(i))) - 1;
        n_samples = floor(segment_length*sample_rate/filters(i));
        
        if mod(n_samples, 2) == 1
            
            n_samples = n_samples - 1;
            
        end
        
        data1_seg = data1( 1:n_samples);        
        Gfilter = (Gfilt(filters(i), filter_width(i), length(data1_seg), sample_rate));
        
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
            
            tap = tukeywin(n_samples, 1);
                                 
            seg1  = real(ifft(Gfilter .* fft(real(seg1).*tap)));
            seg2  = real(ifft(Gfilter .* fft(real(seg2).*tap)));
            
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



%     for i = 1:(length(filters))
%                                       
% %         data1F = butterworthfilt(data1.*tukeywin(length(data1), 0.01),1/sample_rate, ...
% %             filters(i) - filter_width(i)/2, filters(i) + filter_width(i)/2);
% %         data2F = butterworthfilt(data2.*tukeywin(length(data2), 0.01),1/sample_rate, ...
% %             filters(i) - filter_width(i)/2, filters(i) + filter_width(i)/2);
%         tap     = tukeywin(length(data1), 0.01);
%         Gfilter = (Gfilt(filters(i), filter_width(i), length(data1), sample_rate));
% 
%         data1F  = real(ifft(Gfilter .* fft(real(data1).*tap)));
%         data2F  = real(ifft(Gfilter .* fft(real(data2).*tap)));
%         
%         new_sr = min([ sample_rate round(downsample_to*filters(i))]);
% 
%         data1F  = downsample(data1F, round(sample_rate/new_sr));
%         data2F  = downsample(data2F, round(sample_rate/new_sr));
%         
%         %50% overlap
%         n_inc     = floor(length(data1F)/(segment_length*new_sr/filters(i))) - 1;
%         n_samples = floor(segment_length*new_sr/filters(i));
%         
%         for j = 1:n_inc
% 
%             %%%%%%%%%%%%%
%             %pair 1 and 2
%             seg1 = data1F((j-1)*n_samples+1:j*n_samples);
%             seg2 = data2F((j-1)*n_samples+1:j*n_samples);
%                         
%             if integrate
%                
%                 seg1 = cumsum(seg1);
%                 seg2 = cumsum(seg2);
%                 
%             end
%             
%                                     
%             if normalize
%                
%                 seg1 = seg1/rms(seg1);%sign(Z1_seg);%Z1_seg./abs(hilbert(Z1_seg));
%                 seg2 = seg2/rms(seg2);%sign(Z2_seg);%Z2_seg./abs(hilbert(Z2_seg));
%                 
%             end
%             
%             if make_complex
%                
%                 seg1 = hilbert(seg1);
%                 seg2 = hilbert(seg2);
%                 
%             end
%                         
%             if onebit
%                 
%                 seg1 = sign(seg1);
%                 seg2 = sign(seg2);
%                 
%             end
%             
%             cc = corrcoef(seg1,seg2);
%             Call(j) = cc(1,2);
%             
%         end
%                 
%         C(i)       = tanh(median(atanh(real(Call)))) + 1i*median(imag(Call));
%         
%         %get error by bootstrapping
%         C_error(i) = std(bootstrp(bts, @(x) tanh(median(atanh(real(x)))), Call));
%                  
%         Call = [];
%         
%     end
%         
% end

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
