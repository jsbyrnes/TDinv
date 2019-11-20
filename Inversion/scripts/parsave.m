function parsave( modelhist, i, name )
%work around for saving within a parfor loop

    filename = [ name '_chain#' num2str(i) '.mat'];

    save(filename, '-v7.3')
    disp([ 'Saved ' name '_chain#' num2str(i) '.mat'])
    
end

