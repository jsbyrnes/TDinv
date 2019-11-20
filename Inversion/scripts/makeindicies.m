function [ indicies ] = makeindicies( n, m )
%MAKEINDICIES Make a matrix of permutated indicies of size n x m^n
%Gives all the permutation for a non retangular distribution
%For n = 3, gives (1,1,1), (1,1,2), ..., (1,1,m), (1,2,1), ...

    indicies = zeros(n, m^n);

    for i = 1:n

        number = 1;

        for j = 1:m^n

            if mod(j, m^(n - i)) == 0

                if number ~= m

                    number = number + 1;

                else

                    number = 1;

                end

            end

            indicies(i,j) = number;

        end

    end
    
    %I can't figure out why the first and last are backwords
    %oh well circshift
    
    indicies = circshift(indicies, [0,1]);
    

end

