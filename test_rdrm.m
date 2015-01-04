%Test RD RM evolve;
dft = 0;
dfn = 0;
vxlt = Voxel(1, dft, 2.5, 2, 0.05);
vxln = Voxel(1, dfn, 2.5, 2, 0.05);

%time parameters
tot_t = 0.4;
dt = 0.001;
evolving_t = 0 : dt : tot_t;
Nt = length(evolving_t);

%RD field
phase = pi;
tr = 0.006;
rnt = 0.5;

%off resonance
df = 50;

%plot data
Mt = zeros(3, Nt);
Mn = zeros(3, Nt);

vxlt.init_magnet([0; 0; 1]);
vxln.init_magnet([0; 0; 1]);

for it = 1: Nt
  mxm = vxlt.M(1) * (1 - rnt) + vxln.M(1) * rnt;
  mym = vxlt.M(2) * (1 - rnt) + vxln.M(2) * rnt;
  
  vxlt.RD_RMevolve(mxm, mym, phase, df, tr, dt);
  vxln.RD_RMevolve(mxm, mym, phase, df, tr, dt);
  
  Mt(:, it) = vxlt.M;
  Mn(:, it) = vxln.M;
end

%plot
h = figure(31);
subplot(221);
plot(evolving_t, Mt(1, :), 'r', evolving_t, Mn(1, :), 'k');
legend('Mt dw = 5', 'Mn dw = -5');
xlabel('t(s)');
ylabel('Mx');

subplot(222);
plot(evolving_t, Mt(2, :), 'r', evolving_t, Mn(2, :), 'k');
legend('Mt dw = 5', 'Mn dw = -5');
xlabel('t(s)');
ylabel('My');

subplot(223);
plot(evolving_t, Mt(3, :), 'r', evolving_t, Mn(3, :), 'k');
legend('Mt dw = 5', 'Mn dw = -5');
xlabel('t(s)');
ylabel('Mz');

subplot(224);
plot3(Mt(1, :), Mt(2, :), Mt(3,:), 'r', Mn(1, :), Mt(2, :), Mt(3, :), 'k');
legend('Mt dw = 5', 'Mn dw  = -5');
xlabel('Mx');
ylabel('My');
zlabel('Mz');