 function [ solution, s_z, s_vs, s_vpvs, s_fit, tmpbestarray] = ga_drive( HV, Model, ga)
 %GA_drive Performs the GA for HV ratios
 
     %realization

     solution = zeros(length(HV.frequency), ga.maxr);

     %%%% for debugging - save the best rms at each generation for each
     %%%% realization

     for real = 1:ga.maxr

         %initialize the population
         %tmpbestarray = zeros(ga.maxg, 1);
         fprintf('Realization Number %i\n', real);
         
         [members_HV, m_z, m_vs, m_vpvs] = ini_population_HV(ga.population, Model, HV);
         
         %preallocate/clear members and offspring vectors
         for gen = 1:ga.maxg

             fprintf('     Generation Number %i\n', gen);

             %evaluate members
             quality = evaluation(members_HV, HV);

             %%%tmp for testing
             tmpbestarray(gen, 1) = min(quality);

             %create meta data for the offspring
             %for efficiency, I only create them on reproduction
             [o_z, o_vs, o_vpvs] = mate(ga.population, length(Model.z), quality, ga.e_select, m_z, m_vs, m_vpvs);
             
             %disp(num2str( length(unique(o_z(:)) )));
             
             [o_z, o_vs, o_vpvs] = mutate(o_z, o_vs, o_vpvs, ga.mutate, Model);

             %disp(num2str( length(unique(o_z(:)) )));
             
             [members_HV, m_z, m_vs, m_vpvs] = replacement(members_HV, quality, o_z, o_vs, o_vpvs,...
                 ga.e_rep, m_z, m_vs, m_vpvs, HV, Model);

             %disp(num2str(unique(quality)));
             
             %intrarun_results(real, gen) = min(quality);

         end

         semilogy(1:gen, tmpbestarray);
         title('Improvement?')
         
         quality = evaluation(members_HV, HV);

         %bestarray(:, real) = tmpbestarray;

         [~, best] = min(quality);

         solution(:, real) = members_HV(:, best);
         s_z(:, real) = m_z(:, best);
         s_vs(:, real) = m_vs(:, best);
         s_vpvs(:, real) = m_vpvs(:, best);
         s_fit(1, real) = quality(best);

     end
 
 end
 
