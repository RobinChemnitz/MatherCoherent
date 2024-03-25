function [id2m] = mode_array(d_dim, d_shape, d_depth, p_dim, p_shape, p_depth)   
    % Creates a set of selected Fourier modes as the product of two closed
    % balls with respect to a norm and radius specified. 
    % 
    % The input is:
    % d_dim: The dimension of the driving space
    % d_shape: One of the following: 'circle', 'diamond', 'square'. These
    % three shapes correspond to the euclidean norm, the 1-norm and the 
    % infinity-norm.
    % d_depth: The radius of the closed ball for the driving space.
    % p_dim, p_shape, p_depth: The analogous parameters for physical space
    %
    % The output id2m is an N x (d_dim+p_dim) array that contains the
    % information of the N selected Fourier modes. 

    I = -d_depth : d_depth;
    C = cell(1, d_dim);
    for i=1:d_dim
        C{i} = I;
    end
    d_modes = transpose(lattice(C));

    if d_shape == "circle"
        rad = sum(d_modes.^2, 2);
        d_modes = d_modes(rad <= d_depth^2, :);
    end
    if d_shape == "diamond"
        rad = sum(abs(d_modes), 2);
        d_modes = d_modes(rad <= d_depth, :);
    end

    I = -p_depth : p_depth;
    C = cell(1, p_dim);
    for i=1:p_dim
        C{i} = I;
    end
    p_modes = transpose(lattice(C));

    if p_shape == "circle"
        rad = sum(p_modes.^2, 2);
        p_modes = p_modes(rad <= p_depth^2, :);
    end
    if p_shape == "diamond"
        rad = sum(abs(p_modes), 2);
        p_modes = p_modes(rad <= p_depth, :);
    end
    
    d_id2m = kron(d_modes, ones(size(p_modes, 1), 1));
    p_id2m = kron(ones(size(d_modes, 1), 1), p_modes);
    
    id2m = cat(2, d_id2m, p_id2m);
end