%%
%%%%%%%
%This is the driving script for the invesion. It loads and formats data
%from MasterScript_Observation.
%%%%%%%
%%
clear
close all

addpath(genpath('./Inversion/'));

%%%%%%
%define the search
inverse_parameters = define_search;
%%%%%%

%%%%%%
%load the HV ratios and velocities to invert
%The code does not currently handle any components other than single model ZJ0
%[ZR, ZJ0] = load_data(inverse_parameters);

[ZR, ZJ0, model] = syn_data_buriedlayer(inverse_parameters);

%ZR.value     = interp1(ZR.frequency, ZR.value, nf, 'linear', 'extrap');
%ZR.frequency = nf;
ZR.error     = 0.2*ones(size(ZR.error));

%ZJ0.value     = interp1(ZJ0.frequency, ZJ0.value, nf, 'linear', 'extrap');
%ZJ0.frequency = nf;
ZJ0.error     = 0.2*ones(size(ZJ0.value));

%%
p = parpool;
%p.NumWorkers = 1;

%%%%%%%%%%%%%%%%
%Modify this line if 
starting_m = [];

for i = 1:(inverse_parameters.n_startingpoints/p.NumWorkers)

    parfor k = 1:p.NumWorkers
            
        modelhist{:, k} = run_search(ZR, ZJ0, inverse_parameters, k, starting_m);

    end

    if ~exist([inverse_parameters.inversion_name 'inversion1.mat'], 'file')

        savename = [ inverse_parameters.inversion_name 'inversion1.mat'];

    else

        files = dir([inverse_parameters.inversion_name 'inversion*.mat']);
        n = length(files) + 1;

        savename = [ inverse_parameters.inversion_name 'inversion' num2str(n) '.mat'];

    end

    save(savename)
    
    clear modelhist %can be large
    
end
%%%%%%
