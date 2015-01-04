%Simulate tissue contrast under SE and active feedback MRI

%CW parameter
cw_angle = 50;
cw_phi = 0;
tr = 0.002;
rd_fai = 0;

%time parameter
t_tot = 0.1;
dt = 0.001;

t_cw = 1;

%Voxel Group parameter
Nsp = 2000;
T1 = Inf;
T2 = Inf;
df_1 = 20;
df_2 = 20;
df = [df_1, df_2];

dxB_t = 0.7e-6;
dxB_h = 0.8e-6;
dxB = [dxB_t, dxB_h];
d_vessel = 4;
fb_h = 0.02;
fb_t = 0.05;
vblood = [fb_t, fb_h];
D = 1000;
Nvx = 1;
rvxg = 0.5;

%SE parameter
te = 0.1;
N_te = te / (2 * dt);

%vessels parameter
Nvss = 10;
vss_angles = rand(1, 10) * 2 * pi;


sample = Mvv_TVG(Nvx, Nsp, df, T1, T2, Nvss, d_vessel, vblood, vss_angles, D, dxB, rvxg);

h = sample.voxels_1(1).plot_df(0.5, dxB_h);

h1 = sample.test_SEprocess(N_te, t_tot, dt);
filename = 'Tumor Healthy tissue contrast under SE test Rotating Vessels';
print(h1, '-dtiff', filename);

sample = Mvv_TVG(Nvx, Nsp, df, T1, T2, Nvss, d_vessel, vblood, vss_angles, D, dxB, rvxg);
h2 = sample.test_CWprocess(cw_angle, cw_phi, tr, rd_fai, t_cw, dt);
filename = 'Tumor Healthy tissue contrast under CW process rotating vessels';
print(h2, '-dtiff', filename);



