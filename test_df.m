%plot df to test update df functions
%Spin parameters
tic;
Nx = 100;
Ny = 100;
Nspin = Nx * Ny;
df_constant = 33;
v = Voxel(Nspin, df_constant, 1, 4, 0.05);

angle_vessel = 0;
Y_cst = 0.8;

%coordination
x = linspace(0, v.len_space, Nx);
y = linspace(0, v.len_space, Ny);
[X, Y] =meshgrid(x, y);
X = reshape(X, 1, Nspin);
Y = reshape(Y, 1, Nspin);
z1 = ones(1, Nspin) .* v.len_space ./ 4;
z2 = ones(1, Nspin) .* v.len_space;
z3 = ones(1, Nspin) .* v.len_space ./ 2;

pos_1 = [X; Y; z1];
pos_2 = [X; Y; z2];
pos_3 = [X; Y; z3];

v.positions = pos_1;
v.update_bold(angle_vessel, Y_cst);
v.update_dpf();
df_1 = v.df;

v.positions = pos_2;
v.update_bold(angle_vessel, Y_cst);
v.update_dpf();
df_2 = v.df;

v.positions = pos_3;
v.update_bold(angle_vessel, Y_cst);
v.update_dpf();
df_3 = v.df;

toc;

X = reshape(X, Nx, Ny);
Y = reshape(Y, Nx, Ny);
df_1 = reshape(df_1, Nx, Ny);
df_2 = reshape(df_2, Nx, Ny);
df_3 = reshape(df_3, Nx, Ny);
%plot
figure(1);
c1 = contour(X, Y, df_1, 10);
Clabel(c1);
hold on;
v.plot_plane();
title('df at plane of z = len / 4');
xlabel('x/um');
ylabel('y/um');
hold off;

figure(2);
c2 = contour(X, Y, df_2, 10);
Clabel(c2);
hold on;
v.plot_plane();
title('df at plane of z = len');
xlabel('x/um');
ylabel('y/um');
hold off;

figure(3);
C3 = contour(X, Y, df_3, 10);
Clabel(C3);
hold on;
v.plot_plane();
title('df at plane of z = len / 2');
xlabel('x/um');
ylabel('y/um');
hold off;
