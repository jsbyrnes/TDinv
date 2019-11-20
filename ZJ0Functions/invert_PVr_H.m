function [ u_models, fit_history ] = invert_PVr( starting_f, starting_u, uneven_factor, iterations, CZin, f,...
    arclengths, std_search, print_int, burn_in, save_int, chain_ind)
%invert_PVr. Invert the correlation coefficients for rayleigh wave phase
%velocities. This script is for data only taken on the vertical.
%CZ is length(filters) x length(arclengths)

    model_count = 0;

    ffitind = f>=starting_f(1) & f<=starting_f(end);

    nf = length(starting_f);
    nr = length(arclengths);
    
    %interpolate CZ onto starting_f
    for i = 1:length(arclengths)
            
        CZ(i, :) = interp1(f, CZin(i, :), starting_f);
        
    end
        
    xarg = zeros(length(arclengths), nf);
    
    error = rand(1)*0.1;
    
    for i = 1:length(arclengths)
    
        xarg(i,:)            = starting_f*2*pi.*starting_u*arclengths(i);
        bessel_function(i,:) = besselj(0,xarg(i,:))/uneven_factor;
        
    end
            
    starting_fit  = sum((CZ(:) - bessel_function(:)).^2)./(error^2);%chi squared
    fit_history   = zeros(1,iterations);
    current_print = print_int;
    
    for i = 1:iterations
       
        if 100*(i/iterations) > current_print
           
            disp(['Inversion for PVr ' num2str(current_print) '% complete on chain ' num2str(chain_ind) ]);
            current_print = current_print + print_int;
            
        end
                
        index        = randi(length(starting_f) + 1);%randi(length(starting_f) + 1);
        new_u        = starting_u;
        new_error    = error;
        
        if index <= length(starting_f)
        
            p            = normrnd(0, std_search, 1);%perturbation
            new_u(index) = abs(starting_u(index) + p);
            
        else %update the error
           
            new_error = abs(new_error + normrnd(0, 0.01, 1));
            
        end
        
        for k = 1:length(arclengths)
        
            xarg(k, :)            = starting_f*2*pi.*new_u*arclengths(k);
            bessel_function(k, :) = besselj(0,xarg(k, :))/uneven_factor;
            
        end
                
        new_fit                  = sum(((CZ(:) - bessel_function(:)).^2))./(new_error^2);%chi squared
        probability_of_accepting = min([ log(1) ...
            (nr*nf*(log(error/new_error)) - 0.5*(new_fit - starting_fit)) ]);%metropolis law
            
        check = log(rand(1,1));
        
        if probability_of_accepting > check
                    
            starting_fit = new_fit;
            starting_u   = new_u;            
            error        = new_error;
            
        end
            
        %save models
        
        fit_history(i)   = starting_fit;
        error_history(i) = error;
        
        if (i > burn_in) && (mod(i, save_int)==0)
        
            model_count = model_count + 1;
            
            u_models(:, model_count) = starting_u;
            
        end
        
    end

end

