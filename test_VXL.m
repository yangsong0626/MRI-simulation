%test VXL class
%CW paramters
cw_angle = 50;
cw_phi = 0;
tr = 0.002;
rd_fai = 0;

%time parameters
t_tot = 0.05;
dt = 0.0001;

%VXL parameters
Y = 0.5;
T1 = Inf;
T2 = Inf;
df = 0;
D = 1000;
theta = 0;
Nsp = 200;

d_vessel = 4;
vblood = 0.05;%blood volume fraction

%SE parameters
te = 0.04;
N_te = te / dt;

%{
sample_vxl = VXL(Nsp, df, T1, T2, d_vessel, vblood, D, 0);
h1 = sample_vxl.test_CWprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y, false);
filename = 'CW test single VXL 0 degree';
print(h1, '-dtiff', filename);
h2 = sample_vxl.test_SEprocess(N_te, Y, t_tot, dt, false);
filename = 'SE test single VXL 0 degree';
print(h2, '-dtiff', filename);




sample_vxl = VXL(Nsp, df, T1, T2, d_vessel, vblood, D, pi / 2);
h3 = sample_vxl.test_CWprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y, false);
filename = 'CW test single VXL 90 degree';
print(h3, '-dtiff', filename);

h4 = sample_vxl.test_SEprocess(N_te, Y, t_tot, dt, false);
filename = 'SE test single VXL 90 degree';
print(h4, '-dtiff', filename);

h5 = sample_vxl.plot_RW(dt, t_tot);
filename = 'Diffusion Test';
print(h5, '-dtiff', filename);
%}

%{
%Signal Change with d_vessel, vblood is the same
d_vessel = 0: 1: 25;
Nvss = length(d_vessel);
Mzvss = zeros(1, Nvss);
Mxyvss = zeros(1, Nvss);
Mxyse = zeros(1, Nvss);

for ivss = 1: Nvss
  sample_vxl = VXL(Nsp, df, T1, T2, d_vessel(ivss), vblood, D, theta);
  sample_vxl.CWRD_RMprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y, false);
  Mzvss(ivss) = sample_vxl.Mean_Mz;
  Mxyvss(ivss) = sample_vxl.Mxy;
  
  sample_vxl.SE_processI(N_te, Y, t_tot, dt, false);
  Mxyse(ivss) = sample_vxl.Mxy;
end

h = figure;
plot(d_vessel, Mzvss, 'g', d_vessel, Mxyvss, 'b', d_vessel, Mxyse, 'k');
legend('Mz_CW', 'Mxy_CW', 'Mxy_SE');
title('Signal as Vessel diameter changes');
xlabel('Vessel Diameter (um)');
ylabel('M');
filename = 'Signal as vessel diameter changes';
print(h, '-dtiff', filename);
%}

%{
%Signal change as vblood changes.
vblood = 0.01: 0.01: 0.30;
Nvb = length(vblood);
Mzvb = zeros(1, Nvb);
Mxyvb = zeros(1, Nvb);
Mxyse = zeros(1, Nvb);

for ivb = 1: Nvb
  sample_vxl = VXL(Nsp, df, T1, T2, d_vessel, vblood(ivb), D, theta);
  sample_vxl.CWRD_RMprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y, false);
  Mzvb(ivb) = sample_vxl.Mean_Mz;
  Mxyvb(ivb) = sample_vxl.Mxy;
  
  sample_vxl.SE_process(N_te, Y, t_tot, dt, false);
  Mxyse(ivb) = sample_vxl.Mxy;
end
h = figure;
plot(vblood, Mzvb, 'g', vblood, Mxyvb, 'b', vblood, Mxyse, 'k');
legend('Mz_CW', 'Mxy_CW', 'Mxy_SE');
title('Signal as Vessel blood volume fraction changes');
xlabel('Blood Volume Fraction');
ylabel('M');
filename = 'Signal as blood volume fraction changes';
print(h, '-dtiff', filename);
%}

%{
%signal change as Y changes
Y = 0: 0.01: 1;
NY = length(Y);
Mzcw = zeros(1, NY);
Mxycw = zeros(1, NY);
Mxyse = zeros(1, NY);

sample_vxl = VXL(Nsp, df, T1, T2, d_vessel, vblood, D, theta);
for iY = 1: NY
  sample_vxl.CWRD_RMprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y(iY), false);
  
  Mzcw(iY) = sample_vxl.Mean_Mz;
  Mxycw(iY) = sample_vxl.Mxy;
  
  sample_vxl.SE_process(N_te, Y(iY), t_tot, dt, false);
  Mxyse(iY) = sample_vxl.Mxy;
end

h = figure;
plot(Y, Mzcw, 'g', Y, Mxycw, 'b', Y, Mxyse, 'k');
legend('Mz CW', 'Mxy CW', 'Mxy SE');
title('Signal as Oxygen Saturation level changes');
xlabel('Oxygen Saturation Level');
ylabel('M');

filename = 'Signal as Oxygen Saturation Level Y changes';
print(h, '-dtiff', filename);
%}

%Signal Change with theta
theta = 0: 0.1: pi;
Ntheta = length(theta);

Mzcw = zeros(1, Ntheta);
Mxycw = zeros(1, Ntheta);

Mxyse = zeros(1, Ntheta);

for itheta = 1: Ntheta
  sample_vxl = VXL(Nsp, df, T1, T2, d_vessel, vblood, D, theta(itheta));
  
  sample_vxl.CWRD_RMprocess(cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y, false);
  Mzcw(itheta) = sample_vxl.Mean_Mz;
  Mxycw(itheta) = sample_vxl.Mxy;
  
  sample_vxl.SE_process(N_te, Y, t_tot, dt, false);
  Mxyse(itheta) = sample_vxl.Mxy;
end

h = figure;
plot(theta, Mzcw, 'g', theta, Mxycw, 'b', theta, Mxyse, 'k');
legend('Mz CW', 'Mxy CW', 'Mxy SE');
title('Signal as Vessel Angle Changes');
xlabel('Vessel Angle');
ylabel('M');

filename = 'Signal as vessel angle change';
print(h, '-dtiff', filename);

