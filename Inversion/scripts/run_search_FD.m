function run_search( HV, inverse_parameters, name )
%Does the grid search
    
    disp('Setting up the starting points');
    %first, build the set of starting points
    [thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m] = build_startingpoints(inverse_parameters);
    
    %which starting points are alright?
    [model_fit, ~] = evaluate_models(thickness, vs_b, vpvs_b, rho_b, vs_m, vpvs_m, rho_m, HV, inverse_parameters);
    
    [~, sort_ind] = sort(model_fit);
    
    total_kept = 0;
    
    for i = 1:inverse_parameters.neighborhoods
        
        disp([ 'Searching for neighborhood ' num2str(i)]);
        
        %do the ith best starting point
        model_index = sort_ind(i);
        
        disp([ 'Starting fit is ' num2str(model_fit(model_index)) ]);
        
        [nd_fits, nd_HV, nd_thickness, nd_vs_b, nd_vpvs_b, nd_rho_b, ...
            nd_vs_m, nd_vpvs_m, nd_rho_m ] = ...
            neighborhood_search(thickness(:, model_index), vs_b(:, model_index), vpvs_b(:, model_index),...
            rho_b(:, model_index), vs_m(:, model_index), vpvs_m(:, model_index), rho_m(:, model_index), ...
            inverse_parameters, HV, model_fit(model_index));
        
        save_indicies = find(nd_fits < inverse_parameters.threshold);
                
        if ~isempty(save_indicies)
            
            thickness_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_thickness(:, save_indicies);
            vs_b_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_vs_b(:, save_indicies);
            vpvs_b_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_vpvs_b(:, save_indicies);
            rho_b_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_rho_b(:, save_indicies);
            vs_m_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_vs_m(:, save_indicies);
            vpvs_m_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_vpvs_m(:, save_indicies);
            rho_m_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_rho_m(:, save_indicies);
            HV_models_kept(:, total_kept + 1:total_kept + length(save_indicies)) = nd_HV(:, save_indicies);
            model_fits(total_kept + 1:total_kept + length(save_indicies)) = nd_fits(save_indicies);

            total_kept = total_kept + length(save_indicies);
            
        end
        
    end
    
    if total_kept
    
        %%%%%%
        %save the results
        %save the grid of models, the fits, and the predicted ratios
        save(name, 'thickness_kept', 'vs_b_kept', 'vpvs_b_kept', 'rho_b_kept', ...
            'vs_m_kept', 'vpvs_m_kept', 'rho_m_kept','model_fits', 'HV_models_kept', ...
            'inverse_parameters');
    end

end

