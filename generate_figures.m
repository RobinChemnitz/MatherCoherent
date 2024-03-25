% Generate some of the figures seen in the paper.

mkdir('Figures');

% Gyres stream function
v = @(x) [cos(2*pi*x(1, :)).*sin(2*pi*x(2, :)); ...
         -sin(2*pi*x(1, :)).*cos(2*pi*x(2, :))];
psi = @(x) 1/(2*pi) * cos(2*pi*x(1,:)).*cos(2*pi*x(2, :));

res = 40;
I = linspace(0, 1, res);
L = lattice(I, I);
X = reshape(L(1, :), res, res);
Y = reshape(L(2, :), res, res);
stream = psi(L);

fig1 = figure(1);
fig1.Position = [0, 0, 300, 300];
contour(X, Y, reshape(stream, res, res), LineWidth=2);
xticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
yticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
box on
axis square
xlim([0,1]);
ylim([0,1]);
exportgraphics(gcf, 'Figures/Streamfunction.pdf', 'ContentType', 'vector');

res = 17;
I = linspace(0, 1, res);
L = lattice(I, I);
X = reshape(L(1, :), res, res);
Y = reshape(L(2, :), res, res);
dir = v(L);

fig1 = figure(2);
fig1.Position = [0, 0, 300, 300];
quiver(X, Y, reshape(dir(1, :), res, res), reshape(dir(2, :), res, res), LineWidth=1);
xticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
yticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
box on
axis square
xlim([0,1]);
ylim([0,1]);
exportgraphics(gcf, 'Figures/Quiver.pdf', 'ContentType', 'vector');

%%
% Trajectory of oscillation
fig1 = figure(2);
fig1.Position = [0, 0, 300, 300];
alpha = 0.2*[1; sqrt(2)];
delta = 0.15;
func = @(t) delta*[sin(2*pi*t*alpha(1)); cos(2*pi*t*alpha(2))];
T = 10;
t = linspace(0, T, T*100);
offset = func(t);
plot(0.5 + offset(1, :), 0.5 + offset(2, :), LineWidth=1);
xticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
yticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
box on
axis square
xlim([0,1]);
ylim([0,1]);
exportgraphics(gcf, 'Figures/Oscillation.pdf', 'ContentType', 'vector');

%%
mkdir('Figures/Shear');

% Shear vectorfields
v = @(x) 0.5*[sin(2*pi*x(1,:)).*sin(2*pi*x(4, :)); sin(2*pi*x(2,:)).*sin(2*pi*x(3, :))];

theta = [0.25; 0.25];

stream = @(x) 1/2*pi *(-sin(2*pi*theta(1))*cos(2*pi*x(2, :)) + sin(2*pi*theta(2))*cos(2*pi*x(1, :)));



fig = figure(4);
cla();
fig.Position = [0, 0, 300, 300];

res = 40;
I = linspace(0, 1, res);
L = lattice(I, I);
X = reshape(L(1, :), res, res);
Y = reshape(L(2, :), res, res);
Z = reshape(stream(L), res, res);
greymap = 0.7 * ones(2, 3);
contour(X,Y,Z, 5, LineWidth=2);
colormap(greymap);
hold on

res = 17;
I = linspace(0, 1, res);
L = lattice(I, I);
X = reshape(L(1, :), res, res);
Y = reshape(L(2, :), res, res);
Z = reshape(stream(L), res, res);

dir = v([repmat(theta, 1, size(L, 2)); L]);
quiver(X, Y, reshape(dir(1, :), res, res), reshape(dir(2, :), res, res), LineWidth=1);
xticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
yticks([0, 0.2, 0.4, 0.6, 0.8, 1]);
box on
axis square
xlim([0,1]);
ylim([0,1]);
filename = "Figures/Shear/Quiver_025_025.pdf";
exportgraphics(fig, filename, 'ContentType', 'vector');
%close;

