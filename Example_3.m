name = "Shear";
alpha = 0.2*[1; sqrt(3)];
v = @(x) 0.5*[sin(2*pi*x(1,:)).*sin(2*pi*x(4, :)); sin(2*pi*x(2,:)).*sin(2*pi*x(3, :))];
epsilon = 0.03;
v_depth = 1;

G_id2m = mode_array(2, "square", 6, 2, "circle", 8);
G_mode_dic = mode_dictionary(G_id2m);