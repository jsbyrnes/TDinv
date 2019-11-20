function [ kernel ] = get_Vs_kernel( Model, frequency, pert_size)
%get_Vs_kernel Returns the sensitivity of the Model to small Vs perturbations. 
%The pert inc to use is optional, 0.01 km/s if not set. kernel is freq x
%Vs.

    if nargin == 2
        
        pert_size = 0.01;
        
    end

    %get the HV ratios for the current model
    vp = Model.vs.*Model.vpvs;
    starting_HV = get_HV_ratios(Model.z, Model.vs, vp, Model.rho, frequency);
    
    for i = 1:length(Model.vs)
       
        vs_pert = Model.vs;
        
        %perturb each layer to get the derivative
        vs_pert(i) = Model.vs(i) + pert_size;
        pert_HV = get_HV_ratios(Model.z, vs_pert, vp, Model.rho, frequency);
        
        kernel(:, i) = (pert_HV - starting_HV)/pert_size;
        
    end

end

