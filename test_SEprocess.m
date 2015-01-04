df = 0; 
N_spin = 100;
T1 = Inf;
T2 =Inf;
len_vx = 200; 
N_vss = 10; 
d_vss = 20;
sigma = 2;
D = 2.5;
Y = 0.5;

%time properties;
t_tot = 1;
dt = 0.001;
evolving_t = 0: dt: t_tot;
Nt = length(evolving_t);

te = 0.8;
Nte = ceil(te/(2 * dt));

vx = Mvv_VXL(N_spin, df, T1, T2, len_vx, N_vss, d_vss, sigma, D);
[signal_z, signal_xy] = vx.SE_process(Y, dt, Nt, Nte);

%plot
figure(31);
plot(evolving_t, signal_z, 'k', evolving_t, signal_xy, 'r');
legend('Mz', 'Mxy');
