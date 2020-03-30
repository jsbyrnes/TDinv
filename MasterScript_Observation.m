%%%%%%%%%%%%%
%Get data from the IRIS DMC
%Do SPAC and HV measurements for the inversion
%%
clear
clc
close all

addpath('./Functions/');
addpath('./wfTools/');
addpath('./IPGP-sac-matlab-c67a67e');

%%
%%%%%%%%%%%
%Parameters
Parameters = define_parameters( );

%%
%%%%%%%%%%%%%
%get the data
if strcmp(Parameters.file_type, 'SAC')

    dataStruct = load_sac_data(Parameters);
    
elseif strcmp(Parameters.file_type, 'miniseed')
    
    dataStruct = load_miniseed_data(Parameters);
    
elseif strcmp(Parameters.file_type, 'IRIS')
    
    dataStruct = load_data_from_IRIS(Parameters);
    
end
%%
%%%%%%%%
%do SPAC

%make a grid of each station and get the distance between them
list     = unique({dataStruct(:).station});
npairs   = (length(list).^2 - length(list))/2;
station1 = repmat((1:length(list))', [1 length(list)]);
station2 = repmat((1:length(list)), [length(list) 1]);
station1 = tril(station1, -1); station2 = tril(station2, -1);

station_pairs = [ station1(:) station2(:)];
station_pairs = reshape(station_pairs(station_pairs~=0), [ npairs 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do SPAC in the time domain
for n = 1:npairs
    
    [r(n), azi(n)] = distance(dataStruct(station_pairs(n,1)).latitude, dataStruct(station_pairs(n,1)).longitude...
        , dataStruct(station_pairs(n,2)).latitude, dataStruct(station_pairs(n,2)).longitude);
                       
    C = zeros(length(npairs), 3, length(Parameters.freq_range));
    
    for s = 1:length(Parameters.sections)

        T01 = dataStruct(station_pairs(n,1)).T0{s};
        T02 = dataStruct(station_pairs(n,2)).T0{s};
        
        if strcmp(Parameters.time_cull, 'Intersection')

            %cut to the overlapping sections
            [~, i1, i2] = intersect(T01, T02);

        elseif strcmp(Parameters.time_cull, 'EndPoints')

            time_start = max([ T01(1)   T02(1)   ]);
            time_end   = min([ T01(end) T02(end) ]);

            i1 = (T01 > time_start) & (T01 < time_end);
            i2 = (T02 > time_start) & (T02 < time_end);

        end
        
        if length(Parameters.channels) == 3

            data1(:, 1) = double(dataStruct(station_pairs(n,1)).data{1, s});
            data2(:, 1) = double(dataStruct(station_pairs(n,2)).data{1, s});

            %check these!!!
            data1(:, 2) = double(cosd(azi(n))*dataStruct(station_pairs(n,1)).data{2, s} + sind(azi(n))*dataStruct(station_pairs(n,1)).data{3, s});
            data2(:, 2) = double(cosd(azi(n))*dataStruct(station_pairs(n,2)).data{2, s} + sind(azi(n))*dataStruct(station_pairs(n,2)).data{3, s});

            data1(:, 3) = double(-1*sind(azi(n))*dataStruct(station_pairs(n,1)).data{2, s} + cosd(azi(n))*dataStruct(station_pairs(n,1)).data{3, s});
            data2(:, 3) = double(-1*sind(azi(n))*dataStruct(station_pairs(n,2)).data{2, s} + cosd(azi(n))*dataStruct(station_pairs(n,2)).data{3, s});

            if mod(length(data1), 2) == 1

                data1 = data1(1:end - 1, :);
                data2 = data2(1:end - 1, :);

            end

            data1 = data1(i1, :);
            data2 = data2(i2, :);

            for k = 1:length(Parameters.correlations)

                channels_correlate = Parameters.correlations{k};

                %awkward
                if channels_correlate(1) == 'Z'

                    index1 = 1;

                elseif channels_correlate(1) == 'R'

                    index1 = 2;

                elseif channels_correlate(1) == 'T'

                    index1 = 3;

                end

                if channels_correlate(2) == 'Z'

                    index2 = 1;

                elseif channels_correlate(2) == 'R'

                    index2 = 2;

                elseif channels_correlate(2) == 'T'

                    index2 = 3;

                end

                disp([ 'Correlating ' channels_correlate ' for stations ' ...
                    dataStruct(station_pairs(n,1)).station ' & ' dataStruct(station_pairs(n,2)).station ]);

                for i = 1:length(Parameters.segment_length)

                    [ Ctmp(:, i), C_errortmp(:, i) ] = SPAC(data1(:, index1), data2(:, index2),...
                        dataStruct(station_pairs(n,1)).sampleRate, Parameters.segment_length(i), Parameters);

                end
                                
                C(n, k, :)       = C(n, k, :)       + reshape(mean(Ctmp, 2), size(C(n, k, :)));
                C_error(n, k, :) = C_error(n, k, :) + reshape(sqrt(sum(C_errortmp.^2, 2)), size(C(n, k, :)));

            end

        elseif length(Parameters.channels) == 1

            data1(:, 1) = double(dataStruct(station_pairs(n,1)).data{1});
            data2(:, 1) = double(dataStruct(station_pairs(n,2)).data{1});

            data1 = data1(i1, :);
            data2 = data2(i2, :);

            disp([ 'Correlating ZZ for stations ' ...
                dataStruct(station_pairs(n,1)).station ' & ' dataStruct(station_pairs(n,2)).station ]);

            for i = 1:length(Parameters.segment_length)

                [ Ctmp(:, i), C_errortmp(:, i) ] = SPAC(data1, data2,...
                    dataStruct(station_pairs(n,1)).sampleRate, Parameters.segment_length(i), Parameters);

            end

            C(n, 1, :)       = C(n, 1, :)       + reshape(mean(Ctmp, 2), size(C(n, k, :)));
            C_error(n, 1, :) = C_error(n, 1, :) + reshape(sqrt(sum(C_errortmp.^2, 2)), size(C(n, k, :)));

        end

        clear data1 data2 T0
    
    end
    
end

C       = C/length(dataStruct.sections);
C_error = C_error/length(dataStruct.sections);

r = r*111.12*1000;%convert to m

save([ './Observations/SPAC-' Parameters.run_name ], 'Parameters', 'C', 'C_error', 'r', 'azi');%'C_error',
clear C C_error data1 data2 Ctmp C_errortmp index index1 index2 list r azi T0 time_end time_start tind1 tind2

%%
%%%%%%%%
%do ZR

for i = 1:length(dataStruct)
    
    disp(['Measuring ZR for station ' dataStruct(i).station ])
    
    for k = 1:length(Parameters.central_f)
                
        [ R_mean{i, k}, Z_mean{i, k}, T_mean{i, k}, phaseshift{i, k}, ...
            section_length{i, k}, time_start{i, k}, ~, baz_hits{i, k}, CZR{i, k}] = get_ZR(dataStruct(i).data{1}, ...
            dataStruct(i).data{2}, dataStruct(i).data{3}, dataStruct(i).sampleRate, Parameters.central_f(k), ...
            Parameters.halfwidth(k), Parameters.hitlength_cycles, Parameters.baz_step,...
            Parameters.phase_range, Parameters.TR_max, Parameters.max_hits, Parameters.downsample, ...
            Parameters.sections);
        
        time_start{i, k} = time_start{i, k} + dataStruct(i).T0(1);%in datenum from this
                
    end
        
end

save([ './Observations/ZR-' Parameters.run_name ], 'Parameters', 'R_mean', 'Z_mean', 'T_mean', 'time_start', 'baz_hits', ...
    'phaseshift', 'section_length', 'CZR');

%%
%     %now clip
%     T0 = [ dataStruct(station_pairs(n,1)).T0(1) dataStruct(station_pairs(n,2)).T0(1) ];
%     t1 = (((1:dataStruct(station_pairs(n,1)).sampleCount)/dataStruct(station_pairs(n,1)).sampleRate)/(24*60*60) + T0(1));
%     t2 = (((1:dataStruct(station_pairs(n,2)).sampleCount)/dataStruct(station_pairs(n,2)).sampleRate)/(24*60*60) + T0(2));
%         
%     %figure out the overlap times and clip
%     time_start = max( [ max(T0) Parameters.time_window(1) ]);
%     time_end   = min( [ min([t1(end) t2(end)]) Parameters.time_window(2) ]);
%     
%     %pick the start point to align the data on
%     [~, tind1] = min(abs(t1 - time_start));
%     [~, tind2] = min(abs(t2 - time_start));
%     [~, tend1] = min(abs(t1 - time_end));
%     [~, tend2] = min(abs(t2 - time_end));
%     
%     if (tend1 - tind1 + 1) == (tend2 - tind2)
%         
%         tend2 = tend2 - 1;
%         
%     elseif (tend1 - tind1) == (tend2 - tind2 + 1)
%     
%         tend1 = tend1 - 1;
%         
%     end
