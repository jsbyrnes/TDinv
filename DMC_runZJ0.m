%%%%%%%%%%%%%
%Get data from the IRIS DMC
%Do SPAC and HV measurements for the inversion
%%
clear
clc
close all

addpath('./HVfunctions/');
addpath('./ZJ0functions/');
addpath('./wfTools/')
%%
%%%%%%%%%%%
%Parameters

%experiment_name = 'Riverside';
experiment_name = 'Cavola';
%experiment_name = 'TA';

%%%%%%%%%%%%%
%For the data
start_date = '2004-10-10 00:00:00';
end_date   = '2004-10-11 00:00:00';
network    = 'YB';
stations   = { 'CAF5' 'CAI5' }; %'CAE3' ...
    %'CAF5' 'CAF4' 'CAF3' ...
    %'CAG5' 'CAG4' 'CAG3' }; 
vertical   = 'HHZ';
horz_1     = 'HHN';
horz_2     = 'HHE';
location   = '*';

%%%%%%%%%
%for spac
norm_flag      = 1;
complex_flag   = 0;
integrate_flag = 0;
onebit         = 0;
freq_range       = [0.25 5];
df               = 0.1;
filter_width     = 0.1;
segment_length   = 300;

%%%%%%%%%%%%%%
%for ZR ratios
%parameters to control the search for rayleigh waves.
baz_step           = 1; %when testing for a rayleigh wave, the code looks at test back azimuths, searching with this increment.
central_f          = [  1 5 10 ]; %central frequencies to look at. In Hz.
halfwidth          = 0.025;%half width of the Gaussian filter in Hz
phase_range        = 3; %90 degrees +/- this phase shift to for a Rayleigh wave
blip_error_seconds = 0.1; %seconds, for how long can the phase dip out of the range and still be continuous? Sometimes this happens on real rayleigh waves and I don't want this to trigger twice
TR_max             = Inf; %maximum ratio of transverse to radial energy for the test back azimuth. If the transverse is huge the rayleigh wave might not be very good. The results are not very sensitive to this threshold.


%%
%%%%%%%%%%%%%
%get the data
disp('Downloading the data')

if length(stations) == 1

    Z  = irisFetch.Traces(network, stations, location, vertical, start_date, end_date);
    H1 = irisFetch.Traces(network, stations, location, horz_1, start_date, end_date);
    H2 = irisFetch.Traces(network, stations, location, horz_2, start_date, end_date);

else
   
    for k = 1:length(stations)
       
        Z (k) = irisFetch.Traces(network, stations(k), location, vertical, start_date, end_date);
        H1(k) = irisFetch.Traces(network, stations(k), location, horz_1, start_date, end_date);
        H2(k) = irisFetch.Traces(network, stations(k), location, horz_2, start_date, end_date);
        
    end
    
end
    
%now check that you actually got stuff on each channel
z_list  = unique({Z.station});
h1_list = unique({H1.station});
h2_list = unique({H2.station});

kill = zeros(length(z_list), 1);
for k = 1:length(z_list)
   
    if ~any(strcmp(h1_list, z_list(k))) || ~any(strcmp(h2_list, z_list(k)))
    
        kill(k) = 1;
        
    end
        
end

Z      = Z(~kill);%convets to logical
z_list = z_list(~kill);

kill = zeros(length(h1_list), 1);

for k = 1:length(h1_list)
   
    if ~any(strcmp(z_list, h1_list(k))) || ~any(strcmp(h2_list, h1_list(k)))
    
        kill(k) = 1;
        
    end
        
end

H1      = H1(~kill);%convets to logical
h1_list = h1_list(~kill);

kill = zeros(length(h2_list), 1);

for k = 1:length(h2_list)
   
    if ~any(strcmp(z_list, h2_list(k))) || ~any(strcmp(h1_list, h2_list(k)))
    
        kill(k) = 1;
        
    end
        
end

H2      = H2(~kill);%convets to logical
h2_list = h2_list(~kill);

%now check if you need to remove extra samples - not uncommon
%or kill data that is totally the wrong ass length
max_length = 0;

for k = 1:length(z_list)
   
    max_length = max([ max_length length(Z(k).data)]);
    max_length = max([ max_length length(H1(k).data)]);
    max_length = max([ max_length length(H2(k).data)]);
    
end

%if off by one, clip. If off by more than one, kill
kill = zeros(length(z_list), 1);

for k = 1:length(z_list)
   
    if length(Z(k).data) == (max_length)
       
        Z(k).data = Z(k).data(1:(max_length - 1));
        
    elseif length(Z(k).data) < (max_length - 1)
        
        kill(k) = 1;
        
    end
    
    if length(H1(k).data) == (max_length)
       
        H1(k).data = H1(k).data(1:(max_length - 1));
        
    elseif length(H1(k).data) < (max_length - 1)
        
        kill(k) = 1;
        
    end
    
    if length(H2(k).data) == (max_length)
       
        H2(k).data = H2(k).data(1:(max_length - 1));
        
    elseif length(H2(k).data) < (max_length - 1)
        
        kill(k) = 1;
        
    end
    
end

save([ experiment_name '-data']);

%%
%%%%%%%%
%do SPAC

%make a grid of each station and get the distance between them
list     = unique({Z.station});
npairs   = (length(list).^2 - length(list))/2;
station1 = repmat((1:length(list))', [1 length(list)]);
station2 = repmat((1:length(list)), [length(list) 1]);
station1 = tril(station1, -1); station2 = tril(station2, -1);

station_pairs = [ station1(:) station2(:)];
station_pairs = reshape(station_pairs(station_pairs~=0), [ npairs 2]);

%for time domain
f = freq_range(1):df:freq_range(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do SPAC in the time domain
for n = 1:npairs
    
    disp([ 'Time domain SPAC analysis between ' Z(station_pairs(n,1)).station ...
        ' and ' Z(station_pairs(n,2)).station ]);
    
    for i = 1:length(segment_length)

        for j = 1:length(filter_width)

            %disp([ num2str(i) ', ' num2str(j)])
            %[ ftmp, CZ, ~, ~, CZ_error, ~, ~, arclen, sta_azi ] = ...
            %    SPACcycles( D1, [], [], D2, [], [], segment_length(i), filter_width(j), freq_range, norm_flag, -1:0.01:1, 1);%doing 1 for last arg essenitally turns off errors

            [ ftmp, CZ, ~, ~, CZ_error, ~, ~, arclen, sta_azi ] = ...
                SPACcycles( Z(station_pairs(n,1)), [], [], Z(station_pairs(n,2)), [], ...
                [], segment_length(i), filter_width(j), freq_range, ...
                norm_flag, complex_flag, integrate_flag, onebit);%doing 1 for last arg essenitally turns off errors
            
            CZinterp(:,i,j) = interp1(ftmp, CZ, f);

        end

    end

    CZreshape = reshape(CZinterp, length(f), length(segment_length)*length(filter_width));
    
    CZ(n, :) = mean(CZreshape, 2);
    
    %always output, but never changes so you are good
    r(n)   = arclen;
    azi(n) = sta_azi;
    
end 

figure(2)
plot(f, (CZreshape'));
hold on
plot(f, CZ, 'k', 'LineWidth', 3);
xlabel('Hz')
ylabel('CZ')

