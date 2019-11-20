%%%%%%%%%%%%%
%Get data from the IRIS DMC
%Do SPAC and HV measurements for the inversion
%%
clear
clc
close all

addpath('./HVfunctions/');
addpath('./PVrfunctions/');
addpath('./wfTools/')
%%
%%%%%%%%%%%
%Parameters

%%%%%%%%%%%%%
%For the data
start_date    = '2019-06-25 00:00:00';
end_date      = '2019-08-01 00:00:00';
download_inc  = 1;%in days
location      = '*';
network       = 'CI';
station       = { 'LRL' }; %LRL 

stack_width   = 1;%in days

%work in late 2019
%{'ADO','ARV','BAK','BAR','BBR','BC3','CHF','CIA','CWC','DAN','DEC','DGR','DJJ','EDW2','FMP','FUR','GLA','GMR','GRA','GSC','HEC','IKP','IRM','ISA','LRL','MLAC','MPM','MUR','MWC','NEE2','OSI','PASC','PASC','PDM','PHL','PLM','RPV','RRX','SBC','SCI2','SDD','SHO','SLA','SMM','SVD','TIN','TUQ','USC','VES','VOG','VTV'}

vertical   = 'BHZ';
horz_1     = 'BHN';
horz_2     = 'BHE';

%parameters to control the search for rayleigh waves.
baz_step           = 1; %when testing for a rayleigh wave, the code looks at test back azimuths, searching with this increment.
central_f          = [  1 5 10 ]; %central frequencies to look at. In Hz.
halfwidth          = 0.025;%half width of the Gaussian filter in Hz
phase_range        = 3; %90 degrees +/- this phase shift to for a Rayleigh wave
blip_error_seconds = 0.1; %seconds, for how long can the phase dip out of the range and still be continuous? Sometimes this happens on real rayleigh waves and I don't want this to trigger twice
TR_max             = Inf; %maximum ratio of transverse to radial energy for the test back azimuth. If the transverse is huge the rayleigh wave might not be very good. The results are not very sensitive to this threshold.

ZZ = [];
ZN = [];
ZE = [];

%%
%%Loop over the downloading inc 

disp([ 'Starting on ' start_date ' and working till ' end_date ]);
day_vec            = datenum(start_date):download_inc:(datenum(end_date) - download_inc);

baz_array = 0:baz_step:(360-baz_step);

for i = 1:length(day_vec) %assumes that the first day downloaded correctly
    
    disp([ 'Working on ' datestr(day_vec(i), 2) ])
    
    Z  = irisFetch.Traces(network, station, location, vertical, ...
        datestr(day_vec(i), 31), datestr(day_vec(i) + download_inc, 31));
    HN = irisFetch.Traces(network, station, location, horz_1, ...
        datestr(day_vec(i), 31), datestr(day_vec(i) + download_inc, 31));
    HE = irisFetch.Traces(network, station, location, horz_2, ...
        datestr(day_vec(i), 31), datestr(day_vec(i) + download_inc, 31));
    
    if length(Z) == 1 && length(HN) == 1 && length(HE) == 1
    
        %now check that you actually got stuff on each channel
        z_list  = unique({Z.station});
        h1_list = unique({HN.station});
        h2_list = unique({HE.station});

        kill = zeros(length(z_list), 1);
        for k = 1:length(z_list)

            if ~any(strcmp(h1_list, z_list(k))) || ~any(strcmp(h2_list, z_list(k)))

                kill(k) = 1;

            end

        end

        Z      = Z(~kill);%convets to logical
        z_list = z_list(~kill);

        for k = 1:length(h1_list)

            if ~any(strcmp(z_list, h1_list(k))) || ~any(strcmp(h2_list, h1_list(k)))

                kill(k) = 1;

            end

        end

        for k = 1:length(h2_list)

            if ~any(strcmp(z_list, h2_list(k))) || ~any(strcmp(h1_list, h2_list(k)))

                kill(k) = 1;

            end

        end

        %now check if you need to remove extra samples - not uncommon
        %or kill data that is totally the wrong ass length
        max_length = 0;

        for k = 1:length(z_list)

            max_length = max([ max_length length(Z(k).data)]);
            max_length = max([ max_length length(HN(k).data)]);
            max_length = max([ max_length length(HE(k).data)]);

        end

        %if off by one, clip. If off by more than one, kill
        for k = 1:length(z_list)

            if length(Z(k).data) == (max_length)

                Z(k).data = Z(k).data(1:(max_length - 1));

            elseif length(Z(k).data) < (max_length - 1)

                kill(k) = 1;

            end

            if length(HN(k).data) == (max_length)

                HN(k).data = HN(k).data(1:(max_length - 1));

            elseif length(HN(k).data) < (max_length - 1)

                kill(k) = 1;

            end

            if length(HE(k).data) == (max_length)

                HE(k).data = HE(k).data(1:(max_length - 1));

            elseif length(HE(k).data) < (max_length - 1)

                kill(k) = 1;

            end

        end
        
        if mod(length(Z.data), 2) == 1

            Z.data  = Z.data(1:end-1);
            HN.data = HN.data(1:end-1);
            HE.data = HE.data(1:end-1);

        end

        Z.sampleCount  = length(Z.data);
        HN.sampleCount = length(HN.data);
        HE.sampleCount = length(HE.data);
        
    else
        
        kill = 1;
        
    end
    
    if ~any(kill)
               
        for k = 1:length(central_f)
                    
            [ HV_mean{i, k}, ~, ~, ~, ~, R_mean{i, k}, Z_mean{i, k}, T_mean{i, k}, ~, ...
                ~, ~, time_start{i, k}, ~, baz_hits{i, k} ] = get_rayleigh_HV(Z, ...
                HN, HE, central_f(k), halfwidth, (1/central_f(k))*1.5, ...
                blip_error_seconds, baz_array, phase_range, TR_max);
        
            time_start{i, k} = time_start{i, k} + day_vec(i);%in datenum from this
            
        end
            
    elseif any(kill) && i == 1
        
        error('Starting day invalid');
                
    end
    
end

save([network '.' station{:} '.mat' ]);

%%
%%%%%%%%
%do raw correlations
dateaxis   = datenum(start_date):stack_width:(datenum(end_date) - stack_width); 
figure(1)
hold on

for i = 1:length(central_f)

    %HVvec = [HV_mean{:, i}];
    Rvec = [R_mean{:, i}];
    Zvec = [Z_mean{:, i}];
    tvec  = [time_start{:, i}];
    
    for j = 1:length(dateaxis)

        %stacked_result(i, j) = exp(mean(log(HVvec(tvec > dateaxis(j) ...
        %    & tvec < (dateaxis(j) + stack_width)))));

        %data = HVvec(tvec > dateaxis(j) ...
        %    & tvec < (dateaxis(j) + stack_width));
        
        %data = HVvec(tvec > dateaxis(j) ...
        %    & tvec < (dateaxis(j) + stack_width));
        data = Rvec(tvec > dateaxis(j) ...
            & tvec < (dateaxis(j) + stack_width))./...
            Zvec(tvec > dateaxis(j) ...
            & tvec < (dateaxis(j) + stack_width));
        
        if ~isempty(data)
        
            btset = bootstrp(1e3, @(x) exp(mean(log(x))), data);

            stacked_result(i, j) = mean(btset);
            stacked_error(i,j)   = std(btset);
            
        else
           
            stacked_result(i, j) = nan;
            stacked_error(i,j)   = nan;
            
        end
        
    end
    
    errorbar(dateaxis, stacked_result(i, :), 2*stacked_error(i, :), 'o')
    legendlabel{i} = [ num2str(central_f(i)) ' Hz' ];
    
end
datetick('x','mm-dd')
xlim([ 737601 737637 ])
legend(legendlabel, 'Location', 'NorthWest')
plot([datenum('2019-07-04 00:00:00') datenum('2019-07-04 00:00:00')], [ 0 5 ], 'k--')
plot([datenum('2019-07-06 00:00:00') datenum('2019-07-06 00:00:00')], [ 0 5 ], 'k--')
ylim([0.25 3])
ylabel('Rayleigh wave ellipticity')

save([network '.' station{:} '.Ellip.mat' ]);
