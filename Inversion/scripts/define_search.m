function [ inverse_parameters ] = define_search( )
%define_grid_search Enter in the parameters to search over
%For each value, define a number of layers(including the half space) and an
%array of parameters for the value. First cell entry is for first
%parameters, and so on.

    inverse_parameters.n_startingpoints   = 12*2; %number of random models to make within ranges
    inverse_parameters.iter               = 1e4; %number of iterations of the neighborhoods

    %how fine should you make the velocity models when you do the
    %search?
    inverse_parameters.dz           = 1;
    inverse_parameters.depth        = 100;
    inverse_parameters.max_knots    = 10;
    inverse_parameters.collapse     = 10;%remove any layer thinner than this
    
    %Define the range that parameters can take
    inverse_parameters.vs        = [ 200 3000 ];%can't escape this
    inverse_parameters.vpvs      = [ 1.6 3 ];
    inverse_parameters.rho       = [ 1.5 2 ];%used in flag is zero
    inverse_parameters.lock_rho  = 1;%don't solve for rho
    
    inverse_parameters.sigdepths    = 1;
    inverse_parameters.sigvs        = 50;
    inverse_parameters.sigvpvs      = 0.05;
    inverse_parameters.sigrho       = 0.01;%won't use if locked

end
