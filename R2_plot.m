%replicate graph of R2 vs radius under Spin Echo and gradient echo pulse sequence.

%time parameter
t_tot = 0.1;
dt = 0.001;
%SE parameter
te = 0.1;
Nte = te / (2 * dt);

%GRE parameter

%VXL parameter
Nsp = 2000;
df = 0;
T1 = Inf;
T2 = Inf;
f_blood = 0.02;
Diffusion = 1000;
Y = 0.8e-6;

Nvx = 10;
vessel_angles = rand(1, Nvx) * pi * 2;

%BOLD parameter

d_vessel = 2: 2: 200;
Nd = length(d_vessel);
Mxyt = zeros(1, Nd);
Mxyt_2 = zeros(1, Nd);

matlabpool open 2;
parfor id = 1: Nd
  sample = Mvv_VXL(Nsp, df, T1, T2, Nvx, d_vessel(id), f_blood, vessel_angles, Diffusion);
  sample.SE_process(Y, dt, t_tot, Nte);
  Mxyt(id) = sample.Mxy;
end

parfor id = 1: Nd
  sample_2 = Mvv_VXL(Nsp, df, T1, T2, Nvx, d_vessel(id), f_blood, vessel_angles, Diffusion);
  sample_2.FP_process(Y, dt, t_tot);
  Mxyt_2(id) = sample_2.Mxy;
end

matlabpool close;


R2 = -log(Mxyt) / te;
R2_GRE = -log(Mxyt_2) / te;
radius = d_vessel / 2;
h = figure;
semilogx(radius, R2, 'b', radius, R2_GRE, 'k');
xlabel('Radius(um)');
ylabel('-ln(S)/TE (1/s)');
filename = 'R2 test for single vessel model';
print(h, '-dtiff', filename);
