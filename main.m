% Compute and simulate coherent sets of a nonautonomous flow on a
% multidimensional torus that is driven by a quais-periodic rotation. The
% dimension of the driving space is d_dim, the dimension of the physical 
% space is p_dim.
%
% The parameters required are:
% alpha: d_dim x 1 array indicating the quasiperiodic rotation direction.
% v: (d_dim+p_dim) x n -> p_dim x n function handle of a divergence-free 
% vectorfield. The input is an array of n points in (d_dim+p_dim)-dimensional
% space. The output are the directions of the vectorfield in the n points.
% epsilon: The strength of the noise.
% v_depth: The order up to which Fourier-modes of v should be computed. 
% G_mode_dic: A mode dictionary of the Fourier mode used in the galerkin
% approximation of the generator. Either define a custom dictionary or use
% the functions mode_array and mode_dictionary.
% name: The name of the example as a string.

mkdir(name);

%%
% Constructing the augmented generator G
fprintf("Starting computation with %i modes.\n", size(G_mode_dic{1}, 1));
[G] = augmented_generator(alpha, v, epsilon, G_mode_dic, v_depth);

%%
% Computing the smallest magnitute eigenvalues of the generator G. The
% shift by 1 is needed since G is singular.
fprintf("Computing eigenvalues of sparse matrix of size %i x %i.\n", size(G, 1), size(G, 1));
G = G - speye(size(G));
[vec,ev]=eigs(G,min(length(G),100),'sm'); 
G = G + speye(size(G));
ev=diag(ev);
ev = ev + 1;    
disp("Done computing eigenvalues.")

%%
% Plotting the eigenvalues 
evfig = figure(1);
evfig.Position = [0, 0, 400, 300];
scatter(real(ev), imag(ev), 18, "filled");
xlabel('real part');
ylabel('imaginary part');
box on
title('(Partial) spectrum of the augmented generator');
filename = name + "/Spectrum.pdf";
exportgraphics(evfig, filename, 'ContentType', 'vector');
filename = name + "/Spectrum.fig";
savefig(filename);

%%
% Identify the largest eigenvalue with non-zero real part. A threshold of 
% 10^-12 is used to separate zero and non-zero real part.
tmp = ev;
tmp(abs(real(tmp))< 10^-12) = -Inf;
tmp = real(tmp); 
[ev_sorted, ev_order] = sort(tmp, 1, 'descend');

d_ev = ev(ev_order(1));           % The largest eigenvalue
d_vec = vec(:, ev_order(1));      % The corresponding eigenfunction

%%
% Compute and visualize coherent sets. Visualization only works for
% 2D-vectorfields
coherent_sets(alpha, v, epsilon, d_ev, d_vec, G_mode_dic, "positive", 1, name, false);
