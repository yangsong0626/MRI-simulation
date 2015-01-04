%Test RD RM evolve;
vxlt = Voxel(1, 5, 2.5, 2, 0.05);
vxln = Voxel(1, -5, 2.5, 2, 0.05);

%time parameters
tot_t = 1;
dt = 0.001;
evolving_t = 0 : dt : tot_t;
Nt = length(evolving_t);

%RD field
phase = pi;
tr = 0.003;

%off resonance
df = 50;

%steady state pulse sequence
angle = pi / 6;
tss = 0.01;

%plot data
Mt = zeros(3, Nt);
Mn = zeros(3, Nt);

vxlt.init_magnet([0; 0; 1]);
vxln.init_magnet([0; 0; 1]);


for it = 1: Nt
  mxm = vxlt.M(1) + vxln.M(1);
  mym = vxlt.M(2) + vxln.M(2);
  
  vxlt.M = xrot(angle) * vxlt.M;
  vxln.M = xrot(angle) * vxln.M;
  
  vxlt.RD_evolve(mxm, mym, tr, dt);
  vxln.RD_evolve(mxm, mym, tr, dt);
  
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