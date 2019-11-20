function [ tmp ] = plotrfs( RF )
%PLOTRFS Plots everything from the RF struct

    q = length(RF);
    
    figure(16)
    
    zerot = RF.zerot;
    
    tralen = length(RF(1).trace);
    
    if length(unique(zerot)) == 1
        
        t = make_t(zerot(1), tralen - 1, RF(1).dt);
        
        for i = 1:length(RF)
            
            plot(t', RF(i).trace);
%            plot(RF(i).trace);
            hold on
            
        end
        
        xlabel('Seconds')
        
    else
        
        plot(RF(:).trace);
       
        xlabel('Samples');
    end        

end

