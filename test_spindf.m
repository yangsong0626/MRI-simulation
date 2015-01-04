%test spin df evovling
Nspin = 10;
df_constant = 0;
D = 1;
d_vessel = 4;
r_size = 0.05;

v = Voxel(Nspin, df_constant, D, d_vessel, r_size);

angle_vessel = 0;
Y = 0.8;

t_tot = 1; % in s
dt = 0.01;
evovling_t = 0: dt : t_tot;
Nt = length(evovling_t);

df_tot = zeros(Nspin, Nt);

v.init_positions();

for it = 1: Nt
  v.RW();
  v.update_bold(angle_vessel, Y);
  %v.update_dpf();
  df_tot(:, it) = v.df .* dt;
end

%plot
figure(21);
color_list = ['y', 'm', 'c', 'r', 'g', 'b', 'w', 'k', 'y-', 'm-'];
for idf = 1: Nspin
  plot(evovling_t, df_tot(idf, :), color_list(idf));
  hold on;
end
hold off;
xlabel('t(s)');
ylabel('spin off resonance freq evovling');
title('Spin df BOLD AND dipole');