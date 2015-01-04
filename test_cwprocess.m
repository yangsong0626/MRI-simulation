%2voxels test of CW RD
%Voxel parameters
theta = (pi / 18) * (0 : 9);
x = 0 : 9;
Nvx = 10;
Nsp = 2000;
df_vxg1 = 0;
df_vxg2 = 0;

%vessel
d_vessel = 4;
r_size = 0.05;

T1 = Inf;
T2 = Inf;
D = 2.5;

%RD field parameters
tr = 0.002;
RD_phas = 0;

%BOLD model constant
Y_vxg1 = 0.5;
Y_vxg2 = 0.8;

r12 = 0.5;


%time parameters
tot_t = 0.5; %in s
dt = 0.001; %in s
evolving_t = 0: dt : tot_t;
Nt = length(evolving_t);

%CW field parameters
dfcw = 50;
cw_phi = 0;

%SE parameters
te = 0.5;
N_te = ceil(te / (2 * dt));

vxg = TwoVoxelGroup(Nvx, x, Nsp, [df_vxg1, df_vxg2], T1, T2, d_vessel, [r_size, r_size], [Y_vxg1, Y_vxg2], D, theta, r12);

[Mxy1, Mxy2] = vxg.CWRD_RMprocess(dfcw, cw_phi, tr, RD_phas, tot_t, dt, [true, true]);

figure;
plot(evolving_t, Mxy1, 'r', evolving_t, Mxy2, 'k');