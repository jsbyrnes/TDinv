function [ HV_mean, HV_median, Czr_bounded, TV_mean, TR_mean, ...
    R_mean, Z_mean, T_mean, phaseshift_mean, section_length, std_section, ...
    time_start, time_end, baz_hits ] = get_rayleigh_HV_gpu( Z_data, h1_data, h2_data, ...
    central_f, halfwidth, hit_length_seconds, blip_error_seconds, baz_array, phase_range, ...
    TR_max, start_tmp_num )
%get_rayleigh_HV. This is the main driving script.

    HV_mean = zeros(1, 'gpuArray');%initialize but note that they will grow, a lot
    HV_median = zeros(1, 'gpuArray');%initialize but note that they will grow, a lot
    R_mean = zeros(1, 'gpuArray');
    Z_mean = zeros(1, 'gpuArray');
    TV_mean = zeros(1, 'gpuArray');
    T_mean = zeros(1, 'gpuArray');
    phaseshift_mean = zeros(1, 'gpuArray');
    section_length = zeros(1, 'gpuArray');
    std_section = zeros(1, 'gpuArray');
    time_start = zeros(1, 'gpuArray');
    time_end = zeros(1, 'gpuArray');
    baz_hits = zeros(1, 'gpuArray');
    TR_mean = zeros(1, 'gpuArray');
    Czr = zeros(1, 'gpuArray');
    Czr_bounded = zeros(1, 'gpuArray');
        
    HV_count = 0;

    %send the data and needed variables to the gpu
    Z_data.data = gpuArray(complex(single(Z_data.data)));
    h1_data.data = gpuArray(complex(single(h1_data.data)));
    h2_data.data = gpuArray(complex(single(h2_data.data)));
    %     Z_data.sampleRate = gpuArray(Z_data.sampleRate);
    %     h1_data.sampleRate = gpuArray(h1_data.sampleRate);
    %     h2_data.sampleRate = gpuArray(h2_data.sampleRate);

    Z_data.data = Z_data.data - mean(Z_data.data);
    h1_data.data = h1_data.data - mean(h1_data.data);
    h2_data.data = h2_data.data - mean(h2_data.data);

    %I could see this be extremely important for measuring the
    %instentaneous phase
    Z_data.data = detrend(Z_data.data);
    h1_data.data = detrend(h1_data.data);
    h2_data.data = detrend(h2_data.data);

    %divide each trace by the channel's gain
    Z_data.data  = Z_data.data/Z_gain;
    h1_data.data = h1_data.data/N_gain;
    h2_data.data = h2_data.data/E_gain;

    tap = tukeywin(Z_data.sampleCount, 0.9);

    %filter

    Gfilter = gpuArray(Gfilt(central_f, halfwidth, Z_data.sampleCount, (Z_data.sampleRate)));

    %Gfilter = (Gfilt(central_f, halfwidth, Z_data.sampleCount, (Z_data.sampleRate)));

    Z_data.data = real(ifft(Gfilter .* fft(real(Z_data.data).*tap)));
    h1_data.data = real(ifft(Gfilter .* fft(real(h1_data.data).*tap)));
    h2_data.data = real(ifft(Gfilter .* fft(real(h2_data.data).*tap)));

    %only need to do this once
    Z_analytic = hilbert(Z_data.data);
    Z_envelope = abs(Z_analytic);
    Z_phase    = unwrap(angle(Z_analytic))*(180/pi);

    h1_analytic = hilbert(h1_data.data);
    h2_analytic = hilbert(h2_data.data);

    %convert the hit length value to samples
    hit_length = hit_length_seconds*Z_data.sampleRate;
    blip_error = blip_error_seconds*Z_data.sampleRate;

    %get time for the sections, in days
    t = 0:1/Z_data.sampleRate:(Z_data.sampleCount - 1)/Z_data.sampleRate;
    t = t/(24*60*60);

    %rotate into hypothetical R and T
    for i = 1:length(baz_array)

        %rotate into hypothetical R, don't need T
        R_analytic = cosd(baz_array(i))*h1_analytic + sind(baz_array(i))*h2_analytic;
        T_analytic = -1*sind(baz_array(i))*h1_analytic + cosd(baz_array(i))*h2_analytic;

        %get the phase and envelope. Need to unwrap the phase. Radial should be advanced
        %tmp = hilbert(h1_data.data);
        R_envelope = abs(R_analytic);
        R_phase    = unwrap(angle(R_analytic))*(180/pi);

        %tmp = hilbert(h1_data.data);
        T_envelope = abs(T_analytic);
        %T_phase    = unwrap(angle(R_analytic))*(180/pi);

        phase_shift = wrapTo180(R_phase - Z_phase - 90);

        %whereever the phase diffference is within a certain tolerance, get the
        %amplitude ratio
        ind_hit = abs(phase_shift) < phase_range;

        %for each hit, check how long of a continuous segment the phase is
        %consistent
        ind_hit_int = find(ind_hit);
        edges   = diff(ind_hit_int);
        edges = edges > blip_error;
        edges_ind = find(edges);
        section_lengths = diff(edges_ind);
        section_hits = section_lengths >= hit_length;
        sections_use = find(section_hits);

        for q = 1:length(sections_use)

            %section length points to edges_ind
            point = edges_ind(sections_use(q));
            %edges poits to ind_hit_int, add 1
            point = ind_hit_int(point + 1);
            %ind_hit_int points to ind_hit, which matches the data

            length_of_section = section_lengths(sections_use(q)) - 1;

            ind_measure = point:(point + length_of_section);

            TR_meantmp = mean(T_envelope(ind_measure)./R_envelope(ind_measure));

            if TR_meantmp > TR_max

                continue %too much energy on the tangential compared to the radial

            end

            %need to check if there is an overlapping section already saved

            time_start_tmp = start_tmp_num + datenum(t(ind_measure(1)));
            time_end_tmp   = time_start_tmp + datenum(length_of_section/(Z_data.sampleRate*24*60*60));

            start_check = time_start_tmp > time_end;
            end_check   = time_end_tmp < time_start;

            overlapping = find(~(start_check | end_check));

            if ~isempty(overlapping)

                %need to pick the best one. The 'best' one is the max of
                %equation (3) of Stachnik

                %Czr_tmp = sum(imag(R_analytic(ind_measure)).*Z_data.data(ind_measure))/sum(Z_data.data(ind_measure).*Z_data.data(ind_measure));
                Czr_tmp = sum(imag(R_analytic(ind_measure)).*Z_data.data(ind_measure))/sqrt(sum(imag(R_analytic(ind_measure)).*imag(R_analytic(ind_measure))).*sum(Z_data.data(ind_measure).*Z_data.data(ind_measure)));

                Czr_overlapping = Czr_bounded(overlapping);

                %compare the TR ratios of overlapping sections
                if any(Czr_tmp < Czr_overlapping)

                    continue %there is another with a larger correlation

                else

                    %this one is better, save over the previous ones
                    HV_mean(overlapping) = [];
                    HV_median(overlapping) = [];
                    TV_mean(overlapping) = [];
                    TR_mean(overlapping) = [];
                    R_mean(overlapping) = [];
                    Z_mean(overlapping) = [];
                    phaseshift_mean(overlapping) = [];
                    section_length(overlapping) = [];
                    std_section(overlapping) = [];
                    time_start(overlapping) = [];
                    time_end  (overlapping) = [];
                    baz_hits(overlapping) = [];
                    Czr(overlapping) = [];
                    Czr_bounded(overlapping) = [];

                    HV_count = length(HV_mean) + 1;
                    ind_save = HV_count;

                end

            else %no overlapping events

                HV_count = HV_count + 1;
                ind_save = HV_count;

            end

            HV_mean(ind_save) = mean(R_envelope(ind_measure)./Z_envelope(ind_measure));
            HV_median(ind_save) = median(R_envelope(ind_measure)./Z_envelope(ind_measure));
            TV_mean(ind_save) = mean(T_envelope(ind_measure)./Z_envelope(ind_measure));
            TR_mean(ind_save) = mean(T_envelope(ind_measure)./R_envelope(ind_measure));
            R_mean(ind_save) = mean(R_envelope(ind_measure));
            Z_mean(ind_save) = mean(Z_envelope(ind_measure));
            T_mean(ind_save) = mean(T_envelope(ind_measure));
            phaseshift_mean(ind_save) = mean(abs(phase_shift(ind_measure)));
            section_length(ind_save) = length_of_section/Z_data.sampleRate*central_f;
            std_section(ind_save) = std(R_envelope(ind_measure)./Z_envelope(ind_measure));
            time_start(ind_save) = time_start_tmp;
            time_end  (ind_save) = time_end_tmp;
            baz_hits(ind_save) = baz_array(i);
            Czr(ind_save) = sum(imag(R_analytic(ind_measure)).*Z_data.data(ind_measure))/sum(Z_data.data(ind_measure).*Z_data.data(ind_measure));
            Czr_bounded(ind_save) = sum(imag(R_analytic(ind_measure)).*Z_data.data(ind_measure))/sqrt(sum(imag(R_analytic(ind_measure)).*imag(R_analytic(ind_measure))).*sum(Z_data.data(ind_measure).*Z_data.data(ind_measure)));

        end

    end

end

