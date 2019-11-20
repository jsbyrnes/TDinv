function [ t ] = make_t( shift, length , dt)
%MAKE_T Make the time vector for reciever functions
%  The shift and length should be in samples

    t = -(shift*dt):dt:(length-shift)*dt;

end

