function [Fc, vid2m] = fourier_modes_auto(v, dim, depth, progress)
    % Computes the Fourier modes of a function v on a dim-dimensional torus.
    %
    % The input is:
    % v: a function handle f = @(x) [...], of a function T^dim -> C^n of
    % which to compute the Fourier modes.
    % dim: the dimension of the torus.
    % depth: the largst absolute value of any mode in mide_dic. 
    % progress: whether or not to display how many percent of the
    % computation is already completed.
    %
    % The output is:
    % Fc: an n x K matrix, where n is the dimension of the codomain of v
    % and K is the number of Fourier modes with coefficient above a threshold.
    % Each row of Fc contains the Fourier coefficients of the corresponding
    % component of v. Which Fourier modes the individual columns refer to
    % is encoded in vid2m.
    % vid2m: a K x dim integer matrix where each row contains the
    % information of a Fourier mode. 

    err = 0.001; % Threshold below which Fourier coefficients will be 
                 % treated as 0. Enhances sparsity at the cost of errors.
    
    res = 3*depth;

    I = linspace(0,1,res+1); I = I(1:end-1);
    C = cell(1, dim);
    for i=1:dim
        C{i} = I;
    end
    G = lattice(C);
    vals = transpose(v(G)); 

    I = -depth:depth;
    C = cell(1, dim);
    for i=1:dim
        C{i} = I;
    end
    c = transpose(lattice(C));
    
    % Memory saving computation
    n = size(vals, 2);
    N = (2*depth+1)^dim;
    Fc = zeros(n, N);
    batchSize = 100;
    k = 1;
    percent = 0;
    if progress
        disp("Compute Fourier modes:")
    end
    while k * batchSize < N
        ids = (k-1)*batchSize + 1 : k*batchSize;
        basis_func = exp(-2*pi*1i*c(ids, :)*G);
        coeffs = 1/size(G,2) * basis_func * vals;
        Fc(:, ids) = transpose(coeffs);
        if progress && floor((batchSize * k) / N * 100 / 5) * 5 > percent
            percent = floor((batchSize * k) / N * 100 / 5) * 5;
            disp(percent + "%")
        end
        k = k+1;
    end
    ids = (k-1)*batchSize + 1 : N;
    basis_func = exp(-2*pi*1i*c(ids, :)*G);
    coeffs = 1/size(G,2) * basis_func * vals;
    Fc(:, ids) = transpose(coeffs);

    idx = sum(abs(Fc), 1) > err;

    Fc = Fc(:, idx);
    vid2m = c(idx, :);

    fprintf("Found %d non-zero Fouriermodes. Thats %.4f%%. The maximum" + ...
        " depth is %d.\n", size(vid2m, 1), 100*size(vid2m, 1)/N, max(abs(vid2m), [], 'all'));
end