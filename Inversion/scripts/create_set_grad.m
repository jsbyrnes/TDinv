function [ parameters, indicies, pointer] = create_set_grad(base, pert, inc )
%CREAT_SET Creates a matrix of all the parameters you want to use to create
%the model, and the permutation indicies used to go through the matrix

    q = find(pert);
    w = find(inc);

    if q ~= w

        error('Something wrong with the pert and inc arrays - should line up');

    end

    %NaN are ok here
    paramvec = max(pert(q)./inc(w));

    parameters = NaN(length(q), paramvec + 1);

    for i = 1:length(q)

        ind = q(i);

        %This is the only thing that shold be different
        p = base(ind):inc(ind):base(ind) + pert(ind);

        parameters(i, 1:length(p)) = p;

    end

    %make permuting indicies here

    [n m] = size(parameters);

    indicies = makeindicies(n,m);

    %parameter matrix has nans in that don't mean anything
    %find indicies of the nans in parameters and remove those columns
    %from indicies matrix

    for i = 1:n

        x = find(isnan(parameters(i, :)));

        if ~isempty(x)

            for j = length(x):-1:1

                indicies(:, indicies(i, :) == x(j)) = [];

            end

        end

    end
    
    %gives the indicies of the layers that are changing
    pointer = q;
    
    if isempty(indicies)
        
        indicies = 0;
        
    end

    %There you have it - all the necessary permutations
    %Number of created models based on this parameter is the length of indicies

end