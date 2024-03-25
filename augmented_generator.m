function [G] = augmented_generator(alpha, v, epsilon, G_mode_dic, v_depth)
    % Compute a discretitzation of the augmented generator of the Mather
    % semigroup. The input is:
    % alpha: n_d x 1 array indicating the quasiperiodic rotation direction.
    % v: (n_d+n_p) x n -> n_p x n function handle of a divergence-free 
    % vectorfield. 
    % G_mode_dic: A mode dictionary of the Fourier mode used in the galerkin
    % approximation of the generator.
    % v_depth: The order up to which Fourier-modes of v should be computed. 
    %
    % The output is G, a sparse square complex matrix which is a discrete 
    % generator. The size of G is given by the number of Fourier modes in 
    % G_mode_dic.

    Gid2m = G_mode_dic{1};
    Gm2id = G_mode_dic{2};
    
    total_dim = size(Gid2m, 2);
    d_dim = length(alpha);
    p_dim = total_dim - d_dim;

    alpha = reshape(alpha, [d_dim, 1]);
    
    [vF, vid2m] = fourier_modes_auto(v, total_dim, v_depth, true);
    
    N = size(Gid2m, 1);      % Total nuber of modes in G
    M = N;                   % Constant used for memory increments

    diffv = zeros(1,N);
    advi = zeros(1, M);
    advj = zeros(1, M);
    advv = zeros(1, M);

    counter = 1;
    function add_value(a,b,v)
        if v ~= 0
            advi(counter) = a;
            advj(counter) = b;
            advv(counter) = v;
    
            if counter == M
                tmpi = advi;
                tmpj = advj;
                tmpv = advv;
                advi = zeros(1, M+N);
                advj = zeros(1, M+N);
                advv = zeros(1, M+N);
                advi(1:M) = tmpi;
                advj(1:M) = tmpj;
                advv(1:M) = tmpv;
    
                M = M + N;
            end
            counter = counter + 1;
        end
    end

    disp('Construct Advection matrix')
    percent = 0;
    for a = 1:N
        diffv(a) = -0.5*epsilon^2*(2*pi)^2 * sum(Gid2m(a, d_dim+1:end).^2);

        v = Gid2m(a, 1:d_dim) * alpha;
        add_value(a,a,v);

        for k=1:size(vid2m, 1)
            bkey = jsonencode(Gid2m(a, :) - vid2m(k, :));
            if isKey(Gm2id, bkey)
                b = Gm2id(bkey);
                v = Gid2m(b,d_dim+1:end) * vF(:, k);
                add_value(a,b,v);
            end
        end
        
        if floor(a / N * 100 / 5) * 5 > percent
            percent = floor(a / N * 100 / 5) * 5;
            disp(percent + "%");
        end
    end

    D = sparse(1:N, 1:N, diffv);
    advv = -2*pi*1i*advv;
    counter = counter - 1;
    A = sparse(advi(1:counter), advj(1:counter), advv(1:counter));

    G = D + A;
end