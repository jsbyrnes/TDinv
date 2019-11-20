function [ z, vp, vpvs ] = compile_velocity_model( Greens )
%COMPILE_VELOCITY_MODEL Returns the complete velocity structure from the
%Greens struct. Basically, this just moves things around and applies the
%gradient.

    if length(Greens) ~= 1
        
        error('Only supposed to give one Greens struct at a time');
        
    end

    %hardwired. Easy to change though.
    inc = 25;
    
    vp = Greens.vp;
    z = Greens.z;
    vpvs = Greens.vpvs;
    grad = Greens.grad;
        
    if ~all(grad)
        
        for index = 1:length(grad)
            
            if grad(index) && (vp(index) ~= vp(index + 1))
                
                vp_inc = (vp(index + 1) - vp(index))/inc;
                
                %bounds of gradient in depth
                
                %shallow
                b1 = z(index) - grad(index)/2;
                %deeper
                b2 = z(index) + grad(index)/2;
                
                z_inc = (b2 - b1)/inc;
                
                slope_vp = vp(index):vp_inc:vp(index+1) - vp_inc;
                slope_depths = b1:z_inc:b2 - z_inc; % the last depth_inc increment is given by depth(index+1)
                
                if length(slope_vp) ~= length(slope_depths)
                    
                    %shouldn't happen
                    error('slope_vp and slope_depths must be the same length!!');
                    
                end
                
                vp = [vp(1:index - 1) slope_vp vp(index+1:end)];
                z = [z(1:index - 1) slope_depths z(index+1:end)];
                
                q = length(slope_vp);
                
                vpvs = [vpvs(1:index-1) vpvs(index)*ones(1,q) vpvs(index + 1:end)];
                
            end
            
        end
        
    end
    
end

