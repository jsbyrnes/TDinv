function [ P_comp SV_comp] = ZR2PSV(phase, z_comp, r_comp, rp, surfv)
%ZR2PSV Rotate the Z and R traces into ray coordinates
%   surfv should be the vp(1) for a P and vs(1) for an S
%   Make sure that rp is plane wave(s/km). Source component is the 
%   channel that the incoming wave in on, conv component is the channe;
%   that the converted phases are on. Also renormalizes after rotation
%   to incoming component. Assumes all later phases are less than incoming
%   pulse(probably true for synthetic data, want to be true for real data)
%
%   Joseph Byrnes
%   Oct 15th 2012, jbyrnes@uoregon.edu
%

    angle = asin(surfv*rp);

    SV_comp = cos(angle)*r_comp - sin(angle)*z_comp;
    P_comp = sin(angle)*r_comp + cos(angle)*z_comp;
    
    if strcmpi('SV', phase)
    
        amp = max(SV_comp);
    
    elseif strcmpi('P', phase)
        
        amp = max(P_comp);
        
    end
    
    SV_comp = SV_comp/amp;
    P_comp = P_comp/amp;
     
end