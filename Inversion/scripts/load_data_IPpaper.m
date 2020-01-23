function [ HV, PVr, GVr, PVl, GVl ] = load_data_IPpaper( dataname )
%load_data Load data. 
%Possible data are HV ratios (HV), phase velocity of rayleigh waves (PVr),
%group velocity of rayleigh waves (GVr), phase velocity of love waves
%(PVl), and group velocity of low waves (GVl)

    load(dataname, 'HV');
    
    %if you don't have some of these then make the variable empty
    if ~exist('HV');
        
        HV = [];
        
    end

    if ~exist('PVr');
        
        PVr = [];
        
    end

    if ~exist('GVr');
        
        GVr = [];
        
    end

    if ~exist('PVl');
        
        PVl = [];
        
    end

    if ~exist('GVl');
        
        GVl = [];
        
    end

    %check dimensions
    if iscolumn(HV.frequency)
       
        HV.frequency = HV.frequency';
        HV.error     = HV.error';
        HV.value     = HV.value';
        
    end
    
    if iscolumn(PVr.frequency)
       
        PVr.frequency = PVr.frequency';
        PVr.error     = PVr.error';
        PVr.value     = PVr.value';
        
    end
    
end

