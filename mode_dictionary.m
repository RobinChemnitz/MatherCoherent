function [mode_dic] = mode_dictionary(id2m)
    % Creates a container object that encodes a selection of Fourier modes.
    % The input id2m is an N x dim array which contains the information of
    % N selected Fourier modes in dim dimensions.
    % The output mode_dic is a cell array with two entries. The first one
    % is id2m which translates an index to a dim-dimensional Fourier mode.
    % The second entry is a dictionary whose keys are Fourier modes encoded
    % as strings, e.g. '[1,-2,0,5]', and whose values are the corresponding
    % indices in id2m. 

    id2m = cast(id2m, 'double');
    N = size(id2m, 1);
    keys = strings(1, N);

    for k=1:N
       keys(k) = jsonencode(id2m(k, :));
    end
    mode_dic = {id2m, dictionary(keys, 1:N)};
end