function [x_new] = logisticgong2004(x_prev, connectivity_mat)%, a, eps)

eps = 0.8;
a = 1.7;

M_i = sum(connectivity_mat');

x_new = (1-eps)*(1-a*x_prev.^2) + (eps./M_i).*connectivity_mat*(1-a*x_prev.^2);

end
