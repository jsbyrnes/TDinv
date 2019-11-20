clear
close all

dataname = './IrishPark25m/run25mday2_20Hz_complex.mat';
addpath(genpath('./scripts/'));
addpath('mat_disperse');
inversion_name = 'TD_IP';

%%%%%%
%define the search
inverse_parameters = define_search;
%%%%%%

%%%%%%
%load the HV ratios and velocities to invert
%The code does not currently handle the love wave part of this
[HV, ZJ0] = load_IP_data( );
%[HV, ZJ0] = syn_data_buriedlayer(inverse_parameters);

%add noise
% ZJ0.value                = normrnd(ZJ0.value, ZJ0.error);
% ZJ0.value(ZJ0.value > 1) = 1;
% 
% HV.value      = normrnd(HV.value, 0.05);
% HV.error      = 0.05*ones(size(HV.value));

%%%%%%
p = parpool;
%p.NumWorkers = 3;

starting_m = [];

% HV.value     = HV.value(1:2:end);
% HV.error     = HV.error(1:2:end);
% HV.frequency = HV.frequency(1:2:end);
% 
% ZJ0.value = ZJ0.value(1:2:end);
% ZJ0.error = ZJ0.error(1:2:end);
% ZJ0.frequency = ZJ0.frequency(1:2:end);

for i = 1:(inverse_parameters.n_startingpoints/p.NumWorkers)

    parfor k = 1:p.NumWorkers
            
        modelhist{:, k} = run_search(HV, ZJ0, inverse_parameters, k, starting_m);

    end

    if ~exist([inversion_name 'inversion1.mat'], 'file')

        savename = [ inversion_name 'inversion1.mat'];

    else

        files = dir([inversion_name 'inversion*.mat']);
        n = length(files) + 1;

        savename = [ inversion_name 'inversion' num2str(n) '.mat'];

    end

    save(savename)
    
    clear modelhist HVhist ZJ0hist
    
end
%%%%%%
