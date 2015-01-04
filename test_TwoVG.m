%test TwoVoxleGroup class

%CW parameter
cw_angle = 50;
cw_phi = 0;
tr = 0.002;
rd_fai = 0;

%time parameter
t_tot = 0.1;
dt = 0.001;

%Voxel Group parameter
Nvx= 1;
x = 0;
Nsp = 2000;
T1 = Inf;
T2 = Inf;
df_1 = 0;
df_2 = 0;
df = [df_1, df_2];

Y_1 = 0.5;
Y_2 = 0.8;
Y = [Y_1, Y_2];
d_vessel = 4;
vblood = [0.05, 0.05];
D = 2.5;
theta = 0;
rvxg = 0.50;

%SE parameter
te = 0.1;
N_te = te / (2 * dt);

%{
sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
h = sample.test_CWprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, [false, false]);
filename = 'CWtest Two Voxle Group 0 degree r050';
print(h, '-dtiff', filename);

theta = pi / 2;
sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
h = sample.test_CWprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, [false, false]);
filename = 'CWtest Two Voxle Group 90 degree r050';
print(h, '-dtiff', filename);

rvxg = 0.01;
sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
h = sample.test_CWprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, [false, false]);
filename = 'CWtest Two Voxle Group 90 degree r001';
print(h, '-dtiff', filename);

theta = 0;
sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
h = sample.test_CWprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, [false, false]);
filename = 'CWtest Two Voxle Group 0 degree r001';
print(h, '-dtiff', filename);


sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
h1 = sample.test_SEprocess(N_te, t_tot, dt, [false, false]);
filename = 'Two Voxel Group Test SE MRI No rotating with BOLD only 0 degree';
print(h1, '-dtiff', filename);

theta = pi / 2;
sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
h2 = sample.test_SEprocess(N_te, t_tot, dt, [false, false]);
filename = 'Two Voxel Group Test SE MRI No rotating with BOLD only 90 degree';
print(h2, '-dtiff', filename);
%}

%{
%Contrast vs Y (oxygen saturation level)

Y_h = 0.9;
Y_t = 0: 0.01: 0.9;
dY = Y_h - Y_t;
NY = length(Y_t);
contrast_cwz = zeros(1, NY);
contrast_cwxy = zeros(1, NY);
contrast_se = zeros(1, NY);
rotating = [false, false];

for iY = 1: NY
  Y = [Y_h, Y_t(iY)];
  sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
  sample.CWRD_RMprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, rotating);

  contrast_cwz(iY) = abs(sample.Sz_1 - sample.Sz_2);
  
  contrast_cwxy(iY) = abs(sample.Sxy_1 - sample.Sxy_2);
  
  sample.SE_process(N_te, t_tot, dt, rotating);

  contrast_se = abs(sample.Sxy_1 - sample.Sxy_2);
end

h = figure;
plot(dY, contrast_cwz, 'g', dY, contrast_cwxy, 'b', dY, contrast_se, 'k');
legend('CW Mz', 'CW Mxy', 'SE');
title('Contrast as Oxygen Saturation Level (Y) changes');
xlabel('Difference in Oxygen Saturation Level dY');
ylabel('Contrast');

filename = 'Contrast as Oxygen Saturation Level (Y) Changes 0 degree');
print(h, '-dtiff', filename);

theta = pi / 2;
for iY = 1: NY
  Y = [Y_h, Y_t(iY)];
  sample = TwoVoxelGroup(Nvx, x, Nsp, df, T1, T2, d_vessel, vblood, Y, D, theta, rvxg);
  sample.CWRD_RMprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, rotating);

  contrast_cwz(iY) = abs(sample.Sz_1 - sample.Sz_2);
  
  contrast_cwxy(iY) = abs(sample.Sxy_1 - sample.Sxy_2);
  
  sample.SE_process(N_te, t_tot, dt, rotating);

  contrast_se = abs(sample.Sxy_1 - sample.Sxy_2);
end

h1 = figure;
plot(dY, contrast_cwz, 'g', dY, contrast_cwxy, 'b', dY, contrast_se, 'k');
legend('CW Mz', 'CW Mxy', 'SE');
title('Contrast as Oxygen Saturation Level (Y) changes 90 degree');
xlabel('Difference in Oxygen Saturation Level dY');
ylabel('Contrast');

filename = 'Contrast as Oxygen Saturation Level (Y) Changes 90 degree');
print(h1, '-dtiff', filename);
%}

