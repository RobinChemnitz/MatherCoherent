name = "Translated Gyres";
alpha = 0.2*[1; sqrt(2)];
v = @(x) [cos(2*pi*(x(3, :) + x(1, :))).*sin(2*pi*(x(4, :) + x(2, :))); ...
        -sin(2*pi*(x(3, :) + x(1, :))).*cos(2*pi*(x(4, :) + x(2, :)))];
epsilon = 0.03;
v_depth = 1;

% Construct a mode dictionary tailored to the vectorfield
m1 = 2;
m2 = 2;
nrad = 11;

M = 10^6;
G_id2m = zeros(M, 4);
counter = 1;

for n1=-nrad:nrad
    for n2=-nrad:nrad
        if (n1^2 + n2^2) <= nrad^2
            for d1=-m1:m1
                for d2=-m2:m2
                    G_id2m(counter, :) = [n1+d1, n2+d2, n1, n2];
                    counter = counter+1;
                end
            end
        end
    end
end
counter = counter - 1;
G_id2m = G_id2m(1:counter, :);
G_mode_dic = mode_dictionary(G_id2m);