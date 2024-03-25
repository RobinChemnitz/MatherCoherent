name = "Oscillating Gyres";
alpha = 0.2*[1; sqrt(2)];
delta = 0.15;
v = @(x) [cos(2*pi*(x(3, :) + delta*sin(2*pi*x(1, :)))) ...
             .*sin(2*pi*(x(4, :) + delta*cos(2*pi*x(2, :)))); ...
              -sin(2*pi*(x(3, :) + delta*sin(2*pi*x(1, :)))) ...
             .*cos(2*pi*(x(4, :) + delta*cos(2*pi*x(2, :))))];
epsilon = 0.03;
v_depth = 3;

G_id2m = mode_array(2, "square", 6, 2, "circle", 8);
G_mode_dic = mode_dictionary(G_id2m);
