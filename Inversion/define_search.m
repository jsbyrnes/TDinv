function [ inverse_parameters ] = define_search( )
%define_grid_search Enter in the parameters to search over
%For each value, define a number of layers(including the half space) and an
%array of parameters for the value. First cell entry is for first
%parameters, and so on.
    
    inverse_parameters.data_path             = './Observations/';
    inverse_parameters.data_name             = 'Syn1';
    inverse_parameters.inversion_name        = 'Syn1';
    inverse_parameters.radius_threshold      = 2;%m threshold for a "cluster" of radii to stack

    inverse_parameters.n_startingpoints      = 4; %number of random models to make within ranges. Per script, many indentical scripts may be submitted. 
    inverse_parameters.sig                   = 10;%in percent of the range
    inverse_parameters.iter                  = 1e5; %number of iterations of the neighborhoods
    inverse_parameters.writeOutOn            = 1e2;%how often to print to the screen
    inverse_parameters.fit_markers           = [ -1 0 1];%print when a chain fist passes these[ 1 4 7.5 20 50 100 250 ]
    inverse_parameters.interp_style          = 'linear';%string to pass to interp1 for interpolation style
    inverse_parameters.burnin                = 5e4;%kill chain if you haven't converged here - nonconverging ones take FOREVER to finish
    inverse_parameters.saveEach              = 100;%save the model at this iteration; off the update schedule
    inverse_parameters.weights               = [1 1];%ratio of the rates at which vs is updated verses vpvs
    inverse_parameters.birth_style           = 'perturb';%select prior to birth at random locations (but easy to remove)
                                                    %for small changes, but hard to remove. Select perturb for small
                                                    %changes to the model, but hard to delete nodes. 
    inverse_parameters.squeeze_step          = 0; %1 to use the pert that moves nodes together/apart. Don't use. 
    inverse_parameters.delay_HV              = 0; %set error on HV to near infinite until ZJ0 is fit. Don't use. 
    inverse_parameters.enforce_ZJ0           = 0; %don't let the ZJ0 fit exceed 1.5 of the HV fit. Don't use. 
    inverse_parameters.enforce_fasthalfspace = 0; %make the lowest layer the fastest layer. Not working. 
    inverse_parameters.log_depths            = 0; %log2-spacing of depths?
    inverse_parameters.log_ZR                = 1; %log ZR values

    inverse_parameters.fixed_halfspace.flag   = 1;
    inverse_parameters.fixed_halfspace.vs     = 1500;
    inverse_parameters.fixed_halfspace.vpvs   = 2;
    
    %how fine should you make the velocity models when you do the
    %search?
    inverse_parameters.nz                 = 41;%mat_disp is compiled and hardwired to these at the moment. 
    inverse_parameters.depth              = 400;%in m
    inverse_parameters.max_knots          = 20;
    inverse_parameters.min_knots          = 2;
    inverse_parameters.max_vs_curvature   = Inf;%Inf recommended
    inverse_parameters.max_vpvs_curvature = Inf;%Inf recommended
    
    %Define the range that parameters can take
    inverse_parameters.limits.vs   = [ 50 2000 ];%Bounds on Vs
    inverse_parameters.limits.vpvs = [ 1.2 3 ];
    
    %First define a range of fits. When the fit is above that, you use the
    %max, when under the range, you use the min. 

    %hyperparameters
    %Not properly enabled yet. In the works. 
    inverse_parameters.H.sig_style          = 'fixed';%'fixed' requires entering the errors, 'uniform' is same hierchial error for each period. linear is end to end line for errors
    inverse_parameters.H.sigHV              = 0.05;%0.01;%errors
    inverse_parameters.H.sigZJ0             = 0.05;%0.0025;%errors
    inverse_parameters.H.sigHV_limits       = [ 0.01 0.5];%errors
    inverse_parameters.H.sigZJ0_limits      = [ 0.01 0.5];%errors
    
end
