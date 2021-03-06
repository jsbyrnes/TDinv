function dataStruct = load_sac_data(Parameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    files = dir([ Parameters.data_directory '\*.sac']);

    dataStruct.station = [];

    for k = 1:length(files)

        %%%%%%%%%%%
        %How to get statoin latitude and longitude????
        [d, T0, H]  = rdsac([ files(k).folder '\' files(k).name ]);

        if ~any(strcmp({dataStruct.station}, H.KSTNM))

            if isempty([dataStruct.station])

                index = 1;

            else

                index = length({dataStruct.station}) + 1;

            end

        else

            index = find(strcmp({dataStruct.station}, H.KSTNM));

        end

        %extract common fields
        dataStruct(index).sampleRate     = 1/(round(H.DELTA*1e9)/1e9);
        dataStruct(index).station        = H.KSTNM;
        dataStruct(index).latitude       = H.STLA;
        dataStruct(index).longitude      = H.STLO;
        dataStruct(index).T0{1}          = T0 + ((0:(length(d) - 1))...
            /dataStruct(index).sampleRate)/(60*60*24);%assumed common for each channel

        ind = dataStruct(index).T0{1} >= Parameters.time_window(1) &...
            dataStruct(index).T0{1} <= Parameters.time_window(2);
        
        dataStruct(index).sampleCount    = length(d(ind));%assumed common for each channel
        
        T0                      = dataStruct(index).T0{1};
        dataStruct(index).T0{1} = T0(ind);
        chan_ind                = strcmp(H.KCMPNM, Parameters.channels);

        if any(chan_ind)
            
                dataStruct(index).H{chan_ind}           = H;
                dataStruct(index).data{chan_ind, 1}     = Parameters.sign(chan_ind)*double(d(ind));%can be different lengths
                dataStruct(index).azimuth(chan_ind)     = H.CMPAZ;%common fields
                dataStruct(index).inclination(chan_ind) = H.CMPINC;%common fields
                
        else

            disp(['No data loaded for ' H.KSTNM ' on component ' H.KCMPNM ]);

        end

    end
    
end

