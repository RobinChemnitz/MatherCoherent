function [L] = lattice(varargin)
    % Computes a the coordinates of a lattice of arbitrary dimension and
    % gridpoint distribution. Input can be one of the following types:
    % 1) One input variable C, which is a cell array of d 1D-arrays 
    % I_1,...,I_d of lengths n_1, ..., n_d.
    % 2) d 1D-arrays I_1,...,I_d of lengths n_1, ..., n_d.
    % In both cases, the function creates a d-dimensional lattice with 
    % gridpoints I_i in each of the respective dimensions. 
    % The output L is a d x (n_1*...*n_d) matrix containing the coordinates 
    % of each of the (n_1*...*n_d) points.

    input_size = length(varargin);
    if input_size == 1
        C = varargin{1};
    else 
        C = varargin;
    end
    
    d = length(C);
    if d == 1
        L = C{1};
        return;
    end

    dim = zeros(1,d);
    for k=1:d
        dim(k) = length(C{k});
    end
    L = zeros(d, prod(dim));
    
    for i=1:d
        sz = ones(1, d); sz(i) = dim(i);
        vec = reshape(C{i}, sz);
        sz = dim; sz(i) = 1;
        M = repmat(vec, sz);
        L(i, :) = reshape(M(:), 1, prod(dim));
    end
end