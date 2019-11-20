function [ Greens ] = populate_Greens( passedmodel, Greens, name)
%POPULATE Make the synthetics for the Green's structure and save it as a
%giant fle

    q = length(Greens);
    
    %This part is the slowest part of the whole program
    %I tried parfor on my two core machine, it was ~1 second slower
    
    a = exist(['Greens_' name '.mat'], 'file');
    
    if a == 2
        
        load(['Greens_' name '.mat']);
        
    else
                
        len = length(passedmodel.vp);
            
        A = zeros(1, len);
        B = zeros(1, len);
        C = zeros(1, len);
        phi = zeros(1, len);
        theta = zeros(1, len);

        phase = passedmodel.phase;
        dt = passedmodel.dt;
        
        %for i = 1:q
        parfor i = 1:q
            
            [z, vp, vpvs] = compile_velocity_model(Greens(i));
                        
            vs = vp./vpvs;
            
            %for now, everything is isotropic
            
            rho = .33 + .77*vp;
            
            [Greens(i).P, Greens(i).SV, ~ ] = anirec(phase, dt, Greens(i).rp, 0, vp, vs, z, rho, A, B, ...
                C, phi, theta);
            
            if phase == 'P'
                
                [ ~, zerot ] = max(Greens(i).P);
                
            elseif strcmpi(phase,'SV')
                
                [ ~, zerot ] = max(Greens(i).SV);
                
            end
            
            Greens(i).zerot = zerot;
            
        end
        
        save([ 'Greens_' name], 'Greens');
        
    end
    
end

