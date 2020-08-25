function [dataStruct] = load_data_from_IRIS(Parameters)

    dataStruct.station = [];
    
    %get loop for section of time
    sections = Parameters.time_window(1):Parameters.download_sections:(Parameters.time_window(2) - Parameters.download_sections);
    dataStruct.sections = sections;
    
    for k = 1:length(Parameters.station)
        
        for j = 1:length(sections)
        
            [z, h1, h2] = irisFetch_call(Parameters.network{k}, Parameters.station{k}, ...
                Parameters.location, Parameters.channels, datestr(sections(j), 31), ...
                datestr(sections(j) + Parameters.download_sections, 31));

            if isempty(z) || isempty(h1) || isempty(h2)

                
                if j == 1
                    
                    disp(['Data missing for station ' Parameters.station{k} ', network ' Parameters.network{k}]);
                    break%probably nothing at this station, not a perfect fix
                
                else
                    
                    disp(['No data loaded for a section at station ' Parameters.station{k} ', network ' Parameters.network{k}]);
                    continue
                    
                end
                
            else
                
                disp(['Got data for ' datestr(sections(j), 31) ' at station ' Parameters.station{k} ', network ' Parameters.network{k}]);
                
            end

            if ~any(strcmp({dataStruct.station}, z.station))

                if isempty([dataStruct.station])

                    index = 1;

                else

                    index = length({dataStruct.station}) + 1;

                end

            else

                index = find(strcmp({dataStruct.station}, z.station));

            end

            %extract common fields
            dataStruct(index).sampleRate     = round(z.sampleRate);
            dataStruct(index).station        = z.station;
            dataStruct(index).network        = Parameters.network{k};
            dataStruct(index).latitude       = z.latitude;
            dataStruct(index).longitude      = z.longitude;
            dataStruct(index).T0{j}          = sections(j) + ((0:(length(z.data) - 1))...
                /dataStruct(index).sampleRate)/(60*60*24);%assumed common for each channel

            dataStruct(index).sampleCount    = z.sampleCount;%assumed common for each channel

            dataStruct(index).data{1, j}     = Parameters.sign(1)*z.data;%can be different lengths
            dataStruct(index).azimuth(1)     = z.azimuth;%common fields
            dataStruct(index).inclination(1) = z.dip;%common fields

            dataStruct(index).data{2, j}     = Parameters.sign(2)*h1.data;%can be different lengths
            dataStruct(index).azimuth(2)     = h1.azimuth;%common fields
            dataStruct(index).inclination(2) = h1.dip;%common fields

            dataStruct(index).data{3, j}     = Parameters.sign(3)*h2.data;%can be different lengths
            dataStruct(index).azimuth(3)     = h2.azimuth;%common fields
            dataStruct(index).inclination(3) = h2.dip;%common fields
            
        end

    end

end

