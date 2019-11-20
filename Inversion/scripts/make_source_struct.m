function [Source]=make_source_struct(C1_tot, C2_tot, dt, Model)

    % parameters
    %==========================================================================

    sf_norm = 0;          % normalize minimum phase filter amplitude spectrum? (1 = yes, 0 = no)

    use_xcorr = Model.use_xcorr;        % use the cross correlation between C1 and C2 for the spectral smoothing (1 = yes, 0 = no)

    use_log = Model.use_log;          % use the log of the observed spectra for the spectral smoothing (1 = yes, 0 = no)

    %==========================================================================

    [~, nevt] = size(C1_tot);
    
    % init source struct
    Source = struct('npts',0,'df',0,'sps',0,...
        'min_phase_filt',[],'inv_min_phase_filt',[],'spec_est',[],'td_source_est',[]);
    Source = repmat(Source,[nevt 1]);

    % metadata
    npts  = length(C1_tot);
    sps   = 1/dt;
    df = 1/(npts*1/sps);

%   npts = npts*2;    
%     z  = zeros(npts/2, nevt);
%     C1_tot = [C1_tot;z];
%     C2_tot = [C2_tot;z];

    for ii=1:nevt
        
        C1 = C1_tot(:, ii);
        C2 = C2_tot(:, ii);

        % set source metadata fields
        Source(ii).df    = df;
        Source(ii).sps   = sps;
        Source(ii).npts  = npts;
        
        % calculate stacked min-phase source filter
        %==========================================
        % handle zeros in the spectrum
        spec = fft(C1,[],1);
        ipz  = spec == 0;
        spec(ipz) = eps;

        TempFilt = Kolmogorov(C1,'time','freq')./spec;

        TempFilt_i = 1./TempFilt;

        TempFilt   = mean(TempFilt,2);      % using the mean here as the central estimator
        TempFilt_i = mean(TempFilt_i,2);    % using the mean here as the central estimator

        clear spec

        if sf_norm
            % normalize filters
            TempFilt   = TempFilt./abs(TempFilt);
            TempFilt_i = TempFilt_i./abs(TempFilt_i);
        end

        % record filters
        Source(ii).min_phase_filt     = TempFilt;
        Source(ii).inv_min_phase_filt = TempFilt_i;

        % amplitude spectrum

        if use_xcorr % do the xcorr
            C1 = abs( fft(C1,[],1) .* conj( fft(C2,[],1) ) );
            % C1 = sqrt(C1); % correct for the correlation amp spec squaring of the source function
        else
            C1 = abs( fft(C1,[],1) );
        end

        if use_log % do the log
            ipz = C1 == 0;
            if any(any(ipz))
                C1(ipz) = eps;
            end
            C1 = log(C1);
        end

        % calculate estimated source spectrum
        %====================================

        if Model.use_source_smoothing
            
            if Model.smoothing_style == 1 % smoothing by fitting a spline
                
                freq = (0:npts-1).'.*df;
                
                tol = Model.tol; % tolerance (smoothing parameter)
                
                C1 = single(c2saps(freq,double(C1),tol,freq));
                
            elseif Model.smoothing_style == 2 % smoothing with a gaussian
                
                CT = Model.CT;
                
                n    = ceil(npts/2);
                dt   = 1/sps;
                time = [(0:n-1) (n+1:npts)-(npts+1)].*dt;
                
                Gus = exp(-(time).^2./CT^2).';
                
                %save a scaling factor - the gaussian causes a numerical
                %overflow when you take the exponental. Save the max here,
                %normalize the smoothed result, and then multiple by the
                %constant you save.
                
                old_amp = max(C1);
                old_amp = repmat(old_amp, npts, 1);
                
                Gus = fft(Gus);
                Gus = repmat(Gus,1,nevt);
                
                % convolve using fft
                C1 = real( ifft( fft(C1,[],1).*Gus ,[],1));
                
                new_amp = max(C1);
                new_amp = repmat(new_amp, npts, 1);
                
                C1 = C1./new_amp;
                C1 = C1.*old_amp;
                
            end
        
        end

        if use_log % un-do the log
            C1 = exp(C1);
        end

        if use_xcorr
            C1 = sqrt(C1); % correct for the correlation amp spec squaring of the source function
        end

        % combine recorded traces to form specrtal estimate
        Source(ii).spec_est = median(C1,2);

        % estimate the time domain source function
        %=========================================
        temp_spec = Source(ii).spec_est;
        temp_spec = Kolmogorov(temp_spec,'freq','freq'); % create min-phase signal
        temp_spec = temp_spec.*Source(ii).inv_min_phase_filt; % apply inverse phase filter
        temp_spec = real(ifft(temp_spec)); % estimate time series

        Source(ii).td_source_est = temp_spec; % record

    end
    
end

