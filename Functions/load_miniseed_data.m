function dataStruct = load_miniseed_data(Parameters)
%Miniseed data doesn't have metadata so this is **ad hoc**

    files = dir([ Parameters.data_directory '\*.miniseed']);

    dataStruct.station = [];

    for k = 1:length(files)

        %%%%%%%%%%%
        %How to get statoin latitude and longitude????
        [D, ~] = rdmseed([ files(k).folder '\' files(k).name ]);

        s = strsplit(files(k).name, '_');
        s = strsplit(s{1}, '.');%s(1) is network, s(2) is station, s(3) is channel
        
        if ~any(strcmp({dataStruct.station}, s{2}))

            if isempty([dataStruct.station])

                index = 1;

            else

                index = length({dataStruct.station}) + 1;

            end

        else

            index = find(strcmp({dataStruct.station}, s{2}));

        end

        d = cat(1,D.d);
        
        %extract common fields
        dataStruct(index).sampleRate     = D(1).SampleRate;
        dataStruct(index).station        = s{2};
        
        if strcmp(s{2}, 'FX1')

            dataStruct(index).latitude       = 45.15838;
            dataStruct(index).longitude      = -92.935964;
            
        elseif strcmp(s{2}, 'MB1')
            
            %make the lat lon for the next station, offset of 25 m
            dist       = km2deg(25/1000);
            [lat, lon] = reckon(45.15838, -92.935964, dist, 270);
           
            dataStruct(index).latitude       = lat;
            dataStruct(index).longitude      = lon;
            
        end
        
        dataStruct(index).sampleCount    = length(d);%assumed common for each channel

        chan_ind = strcmp(s{3}, Parameters.channels);

        if any(chan_ind)

            %dataStruct(index).H(chan_ind)           = H;
            dataStruct(index).data{chan_ind}        = d;%can be different lengths
            
            if strcmp(Parameters.channels(chan_ind), 'HHZ')
            
                dataStruct(index).azimuth(chan_ind)     = 0;%common fields
                dataStruct(index).inclination(chan_ind) = 0;%common fields
                
            elseif strcmp(Parameters.channels(chan_ind), 'HHX')
                
                dataStruct(index).azimuth(chan_ind)     = 0;%common fields
                dataStruct(index).inclination(chan_ind) = 90;%common fields
                
            elseif strcmp(Parameters.channels(chan_ind), 'HHY')
                
                dataStruct(index).azimuth(chan_ind)     = 90;%common fields
                dataStruct(index).inclination(chan_ind) = 90;%common fields
                
            end

            dataStruct(index).T0{chan_ind}  = cat(1,D.t);%assumed common for each channel
            
        else

            disp(['No data loaded for ' s{2} ' on component ' s{3} ]);

        end

    end
    
    if length(Parameters.channels) == 1
    
        t1 = dataStruct(1).T0{1};
        t2 = dataStruct(2).T0{1};
            
        d1 = dataStruct(1).data{1};
        d2 = dataStruct(2).data{1};
        
        t2 = t2( (t2 >= Parameters.time_window(1))...
            & (t2 <= Parameters.time_window(2)) );

        %cut to the overlapping sections
        [~, ia, ib] = intersect(t1,t2);
        dataStruct(1).data{1}     = d1(ia);
        dataStruct(2).data{1}     = d2(ib);
        
        dataStruct(1).sampleCount = length(dataStruct(1).data{1});
        dataStruct(2).sampleCount = length(dataStruct(2).data{1});
        
        dataStruct(1).T0 = t1(ia);
        dataStruct(2).T0 = t2(ib);
        
    elseif length(Parameters.channels) == 3
    
        for k = 1:length(dataStruct)

            tZ = dataStruct(k).T0{1};
            t1 = dataStruct(k).T0{2};
            t2 = dataStruct(k).T0{3};

            dZ = dataStruct(k).data{1};
            d1 = dataStruct(k).data{2};
            d2 = dataStruct(k).data{3};

            [~, ia, ib]     = intersect(t1,t2);
            d1              = d1(ia);
            d2              = d2(ib);
            t1              = t1(ia);
            t2              = t2(ib);

            [~, ia, ib]     = intersect(tZ,t1);
            dZ              = dZ(ia);
            d1              = d1(ib);
            tZ              = tZ(ia);

            [~, ia, ib]     = intersect(tZ,t2);
            dZ              = dZ(ia);
            d2              = d2(ib);

            dataStruct(k).data{1} = dZ;
            dataStruct(k).data{2} = d1;
            dataStruct(k).data{3} = d2;

            dataStruct(k).sampleCount = length(dZ);
            dataStruct(k).T0{1}       = tZ;

        end

    end
    
end

%     [D, ~] = rdmseed(data1);
%     samplerate = D(1).SampleRate;
%     datatmp = cat(1,D.d);
%     tZ = cat(1,D.t);
%     Z_data.data = datatmp;
%     Z_data.sampleRate = samplerate;
%     Z_data.sampleCount = length(datatmp);
%     Z_data.latitude  = 45.15838;
%     Z_data.longitude = -92.935964;
% 
%     [D, ~] = rdmseed(data2);
%     samplerate = D(1).SampleRate;
%     datatmp = cat(1,D.d);
%     t1 = cat(1,D.t);
%     h1_data.data = datatmp;
%     h1_data.sampleRate = samplerate;
%     h1_data.sampleCount = length(datatmp);
%     h1_data.latitude  = 45.15838;
%     h1_data.longitude = -92.935964;
% 
%     [D, ~] = rdmseed(data3);
%     samplerate = D(1).SampleRate;
%     datatmp = cat(1,D.d);
%     t2 = cat(1,D.t);
%     h2_data.data = datatmp;
%     h2_data.sampleRate = samplerate;
%     h2_data.sampleCount = length(datatmp);
%     h2_data.latitude  = 45.15838;
%     h2_data.longitude = -91.93596;


%     %cut to the overlapping sections
% 
%     t_start = max([ tZ(1) t1(1) t2(1) ]);
%     t_end   = min([ tZ(end) t1(end) t2(end) ]);
% 
%     [~, ts] = min(abs(t_start - tZ)); 
%     [~, te] = min(abs(t_end - tZ)); 
%     tZn = tZ(ts:te);
%     Z_data.data = Z_data.data(ts:te);
%     [~, ts] = min(abs(t_start - t1)); 
%     [~, te] = min(abs(t_end - t1)); 
%     t1n = t1(ts:te);
%     [~, ts] = min(abs(t_start - t2)); 
%     [~, te] = min(abs(t_end - t2)); 
%     t2n = t2(ts:te);
%     % 
%     % % figure(1)
%     % % hold on
%     % % plot(t1, h1_data.data)
%     % % figure(2)
%     % % hold on
%     % % plot(t2, h2_data.data)
%     %now interpolate, there was a clock issue
%     h1_data.data = interp1(t1, h1_data.data, tZn);
%     h2_data.data = interp1(t2, h2_data.data, tZn);
% 
%     tZn          = tZn(100*10*60:end);
%     Z_data.data  = Z_data.data(100*10*60:end);
%     h1_data.data = h1_data.data(100*10*60:end);
%     h2_data.data = h2_data.data(100*10*60:end);
% 
%     % figure(1)
%     % hold on
%     % plot(tZn, h1_data.data)
%     % figure(2)
%     % hold on
%     % plot(tZn, h2_data.data)
%     % 
%     % figure(3)
%     % plot(tZn, Z_data.data - medfilt1(Z_data.data, 5000))
%     % hold on
%     % plot(tZn, h1_data.data - medfilt1(h1_data.data, 5000))
% 
