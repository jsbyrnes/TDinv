function [ HV, ZJ0 ] = load_IP_data( )
%load_data Load data. 
%Possible data are HV ratios (HV), phase velocity of rayleigh waves (PVr),
%group velocity of rayleigh waves (GVr), phase velocity of love waves
%(PVl), and group velocity of low waves (GVl)

    load('./IrishPark25m/newIRresult.mat', 'HV');
            
    load('./IrishPark25m/run25mday2_20Hz_0.5Hz.mat');
    
    ZJ0.r         = 25;
    ZJ0.value     = real(CZ_TD);
    ZJ0.error     = CZ_TD_error;
    ZJ0.frequency = f';
    
    ZJ0.value     = ZJ0.value(1:2:end);
    ZJ0.frequency = ZJ0.frequency(1:2:end);
    ZJ0.error     = ZJ0.error(1:2:end);
    
    HV.value     = interp1(HV.frequency, HV.value, 0.5:1:19.5, 'linear');
    HV.frequency = 0.5:1:19.5;
    HV.error     = 0.12*ones(size(HV.value));
    
    %check dimensions
    if iscolumn(HV.frequency)
       
        %HV.frequency = HV.frequency';
        HV.error     = HV.error';
        HV.value     = HV.value';
        
    end
    
    if iscolumn(ZJ0.frequency)
       
        ZJ0.frequency = ZJ0.frequency';
        ZJ0.error     = ZJ0.error';
        ZJ0.value     = ZJ0.value';
        
    end
    
end

