function [ P_comp, SV_comp, SH_comp ] = anirec( phase, dt, slow, baz, vp, vs, z, rho, A, B, C, phi, theta, rotate)
%ANIREC Runs the anirec program, and returns a normalized three component green's
%function in ray coordinates. Length of returned trace is fixed to 100 seconds. Phase is case
%insensitive. dt in seconds, slow in sec/km, baz in degree clockwise from
%north. Both velocities in km/s, z is depth to base of layer. Last layers
%is infinite half space. 

    phase = upper(phase);

    if ~strcmpi(phase, 'P') && ~strcmpi(phase, 'SV')
        
        error('Phase should be P or SV');
        
    end
    
    vp = vp*1000;
    vs = vs*1000;
    z = z*1000;
    rho = rho*1000;
    
    cc = 1/slow;
    
    if strcmpi(phase, 'P')        

        [z_comp, r_comp, SH_comp] = anirec_gateway(theta, phi, z, vp, ...
            A, B, vs, C, rho, length(vp) - 1, cc, baz, dt, 1);
        
        if rotate
        
            [P_comp, SV_comp] = ZR2PSV('P', z_comp, r_comp, slow, vp(1)/1000);
            
        elseif ~rotate
            
            P_comp = z_comp;
            SV_comp = r_comp;
            
        end
        
    elseif strcmpi(phase, 'SV')
        
        [z_comp, r_comp, SH_comp] = anirec_gateway(theta, phi, z, vp, ...
            A, B, vs, C, rho, length(vp) - 1, cc, baz, dt, 2);
        
        %z_comp = -1*z_comp;
        
        if max(abs(r_comp)) > max(r_comp)
            
            r_comp = r_comp*-1;
            z_comp = z_comp*-1;
            
        end
        
        if rotate
        
            [P_comp, SV_comp] = ZR2PSV('SV', z_comp, r_comp, slow, vs(1)/1000);
        
        elseif ~rotate
            
            P_comp = z_comp;
            SV_comp = r_comp;
            
        end
                
    end
    
end

