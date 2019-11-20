function [c,u,e] = run_disper80(T,theModel,RorL, maxdepth)
% CALL_DISPER80: Pirated from brandon's code, modified to work with
% standard out. RorL must be 4 for Rayleigh waves or 3 for Love waves

    [ n_layers, ~ ] = size(theModel);

    modelstr = [ ];

    theModel = theModel';
    theModel = theModel(:);
    
    for i = 1:length(theModel)

        modelstr = [ modelstr ' ' num2str(theModel(i), 6) ];

    end    

    Tstr = [];

    for i = 1:length(T)

        Tstr = [ Tstr  ' ' num2str(T(i)) ];

    end

    np = length(T);

    commandstr = [ './DISPER80_src/disper80 ' num2str(n_layers) ' ' modelstr ' ' num2str(np) ' ' ...
        num2str(maxdepth/100) ' ' num2str(RorL) ' ' Tstr ];

    % Call DISPER80--syntax supresses DISPER80 output
    [~,out] = system(commandstr);

    out = strsplit(out);
    out = out(~cellfun('isempty', out));
    
    try
    
    out = cellfun(@str2num, out);
    
    catch
        
        disp('Nonuniform output');
        c = 1e9*ones(size(T));
        u = 1e9*ones(size(T));
        e = 1e9*ones(size(T));
        return
        
    end
        
    if length(out) == 8*length(T)
        
        out = reshape(out, [ 8 length(T) ]);
        c   = out(3,:);
        u   = out(4,:);
        e   = abs(out(6,:).*out(1,:));
        
    else%sometimes periods get skipped, interpolate to the right grid
        
        if mod(length(out),8) == 0
        
            out = reshape(out, [ 8 length(out)/8 ]);
            
        else
            
            %disp('Length of out not right')
            c = 1e9*ones(size(T));
            u = 1e9*ones(size(T));
            e = 1e9*ones(size(T));
            return
            
        end
        
        Tout = out(2,:);
        
        if length(Tout) < 3 || length(unique(Tout)) ~= length(Tout)
            
            %disp('Length of Tout not right')
            c = 1e9*ones(size(T));
            u = 1e9*ones(size(T));
            e = 1e9*ones(size(T));
            return

        end
            
        [Tout, ind] = sort(Tout);
        
        c    = interp1(Tout, out(3,ind), T, 'linear', 'extrap');
        u    = interp1(Tout, out(4,ind), T, 'linear', 'extrap');
        e    = interp1(Tout, abs(out(6,ind).*out(1,ind)), T, 'linear', 'extrap');
        
    end
        
    %remove spikes, these are due to issues in disper80/interp.f
    c = medfilt1(c, 3);
    u = medfilt1(u, 3);
    e = medfilt1(e, 3);
    
end
