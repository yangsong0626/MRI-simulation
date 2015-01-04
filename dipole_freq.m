function DF_freq = dipole_freq(spin_pos, N_spin, dipole_pos, N_dipole, r_dipole, dw_rms)
 %rms angular frequency shift at particle surface in rad/s
  gyro = 267.513e6; %gyromagnetic ratio of spin_pos, in rad/(s*T)
  C = sqrt(5/4);%Constant
  DF_freq = zeros(1, N_spin);
  
  for id = 1: N_dipole
    distance = sqrt((spin_pos(1, :) - dipole_pos(1, id)) .^ 2 + (spin_pos(2, :) - dipole_pos(2, id)) .^ 2 + (spin_pos(3, :) -dipole_pos(3, id)) .^ 2);
    
    for isp = 1 : N_spin
      if distance(isp) < r_dipole
        distance(isp) = r_dipole;
      end
    end
    
    cos_theta = (spin_pos(3, :) - dipole_pos(3, id)) ./ distance;
    freq_temp = 2 .* pi .* C .* r_dipole .* dw_rms .* r_dipole .* (3 .* cos_theta .^ 2 - 1) ./  distance;
    DF_freq = DF_freq + freq_temp;
  end
  
  
    