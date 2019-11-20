function [ Greens ] = populate_Greens( Model, Greens, name)
%POPULATE Make the synthetics for the Green's structure and save it as a
%giant fle

    q = length(Greens);
    
    %This part is the slowest part of the whole program
    %I tried parfor on my two core machine, it was ~1 second slower
    
    a = exist(['Greens_' name '.mat'], 'file');
    
    if a == 2
        
        load(['Greens_' name '.mat']);
        
    else
        
        for i = 1:q
            
            [z, vp, vpvs] = compile_velocity_model(Greens(i));
            
            poisson = vpvs2poisson(vpvs);
            
            writeicmod_3(vp, z, poisson, name, 2);
            
            !chmod a+x icmodinput.csh
            !./icmodinput.csh
            
            [Z,R] = drive_respknt(name, Model.phase, Model.dt, Model.len, Greens(i).rp, 'f' , 'y');
            
            %put in P-SV(returns by respknt as Z-R)
            
            if Model.phase == 'P'
                
                [P, SV] = ZR2PSV(Z, R, Greens(i).rp, vp(1));
                
            elseif Model.phase == 'S'
                
                [SV, P] = ZR2PSV(R, Z, Greens(i).rp, vp(1)/vpvs(1));
                
            end
            
            %return these traces, which are essentially green's functions
            
            %normalize to direct arrival on the P channel for a P, SV for
            %an S
            
            if Model.phase == 'P'
                
                [amp zerot] = max(P);
                
            elseif Model.phase == 'S'
                
                [amp zerot] = max(SV);
                
            end
            
            Greens(i).P = P/amp;
            Greens(i).SV = SV/amp;
            Greens(i).zerot = zerot;
            
        end
        
        save([ 'Greens_' name], 'Greens');
        
    end
    
end

