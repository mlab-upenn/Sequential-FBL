vec_x=[5 1.1 9];
C=30;
% x_min = (-1/C)*log(sum(exp(-C*vec_x)));
x_min = (-1/C)*log(sum(exp(-C*vec_x)));
log(sum(exp(-C*vec_x)))