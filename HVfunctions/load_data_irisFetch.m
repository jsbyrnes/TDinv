function [ Z_data, h1_data, h2_data ] = load_data_irisFetch( sta_network, sta_name, sta_loc, ...
 channel_list, start_tmp, end_tmp, decimate_factor)
%load_data_irisFetch 

    pause(1)
    Z_data = irisFetch.Traces(sta_network, sta_name, sta_loc, channel_list(1, :), start_tmp, end_tmp);
    pause(1)
    h1_data = irisFetch.Traces(sta_network, sta_name, sta_loc, channel_list(2, :), start_tmp, end_tmp);
    pause(1)
    h2_data = irisFetch.Traces(sta_network, sta_name, sta_loc, channel_list(3, :), start_tmp, end_tmp);

    %need to do checks here to see if you got anything
    if isempty(Z_data) || isempty(h1_data) || isempty(h2_data)

        Z_data  = [];
        h1_data = [];
        h2_data = [];

        return

    end

    %For gaps in the data. Fill it in with zeros,
    %for each channel individually.

    n_z = length(Z_data);

    for q = 1:n_z - 1

        gap = round((Z_data(q).endTime - Z_data(q + 1).startTime)*24*60*60/Z_data(1).sampleRate);%convert to second, then samples

        Z_data(1).data = [ Z_data(1).data; zeros(gap, 1); Z_data(q + 1).data ];

    end

    Z_data = Z_data(1);
    Z_data.sampleCount = length(Z_data.data);

    n_h1 = length(h1_data);

    for q = 1:n_h1 - 1

        gap = round((h1_data(q).endTime - h1_data(q + 1).startTime)*24*60*60/h1_data(1).sampleRate);%convert to second, then samples

        h1_data(1).data = [ h1_data(1).data; zeros(gap, 1); h1_data(q + 1).data ];

    end

    h1_data = h1_data(1);
    h1_data.sampleCount = length(h1_data.data);

    n_h2 = length(h2_data);

    for q = 1:n_h2 - 1

        gap = round((h2_data(q).endTime - h2_data(q + 1).startTime)*24*60*60/h2_data(1).sampleRate);%convert to second, then samples

        h2_data(1).data = [ h2_data(1).data; zeros(gap, 1); h2_data(q + 1).data ];

    end

    h2_data = h2_data(1);
    h2_data.sampleCount = length(h2_data.data);


    %sometimes the traces are slightly different sizes for some reason. The
    %end is different between traces.
    if (Z_data.endTime - h1_data.endTime) > 0

        cut = abs(round((Z_data.endTime - h1_data.endTime)*24*60*60/Z_data.sampleRate));
        Z_data.data = Z_data.data(1:end - cut);

        Z_data.endTime = h1_data.endTime;
        Z_data.sampleCount = length(Z_data.data);

    end

    if (Z_data.endTime - h2_data.endTime) > 0

        cut = abs(round((Z_data.endTime - h2_data.endTime)*24*60*60/Z_data.sampleRate));
        Z_data.data = Z_data.data(1:end - cut);
        Z_data.sampleCount = length(Z_data.data);
        Z_data.endTime = h2_data.endTime;

    end

    if (Z_data.endTime - h1_data.endTime) < 0

        cut = abs(round((Z_data.endTime - h1_data.endTime)*24*60*60/Z_data.sampleRate));
        h1_data.data = h1_data.data(1:end - cut);

        h1_data.endTime = Z_data.endTime;
        h1_data.sampleCount = length(h1_data.data);

    end

    if (Z_data.endTime - h2_data.endTime) < 0

        cut = abs(round((Z_data.endTime - h2_data.endTime)*24*60*60/Z_data.sampleRate));
        h2_data.data = h2_data.data(1:end - cut);

        h2_data.endTime = Z_data.endTime;
        h2_data.sampleCount = length(h2_data.data);

    end

    if mod(Z_data.sampleCount, 2) ~= 0

        Z_data.data = Z_data.data(1:end - 1);
        Z_data.sampleCount = length(Z_data.data);

    end

    if mod(h1_data.sampleCount, 2) ~= 0

        h1_data.data = h1_data.data(1:end - 1);
        h1_data.sampleCount = length(h1_data.data);

    end

    if mod(h2_data.sampleCount, 2) ~= 0

        h2_data.data = h2_data.data(1:end - 1);
        h2_data.sampleCount = length(h2_data.data);

    end

    if (Z_data.sampleCount ~= h1_data.sampleCount) || (Z_data.sampleCount ~= h2_data.sampleCount)

        Z_data = [];
        h1_data = [];
        h2_data = [];

        return

    end

    %decimate the traces to a factor what I'll filter to
    %Does not work. 
    if decimate_factor > 0

        %decimate to the new sampling rate
        tmp_data = resam([ Z_data.sampleRate h1_data.sampleRate h2_data.sampleRate] , 80, ...
            [ Z_data.data; h1_data.data; h2_data.data]);
        Z_data.data = tmp_data(:, 1);
        h1_data.data = tmp_data(:, 2);
        h2_data.data = tmp_data(:, 3);

    end

    if mod(length(Z_data.data),2)==1

        Z_data.data = Z_data.data(1:end-1);
        h1_data.data = h1_data.data(1:end-1);
        h2_data.data = h2_data.data(1:end-1);

    end

end

