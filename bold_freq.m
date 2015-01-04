function BOLD_freq = bold_freq(num_spin, d_vessel, theta, spin_pos, Y, space_center)
%generate BOLD field frequency for every spin
d_chi = 4 * pi * 0.277 * 10 ^ -6;
gyro = 42.6e6;
B0 = 7;
Constant = B0 * gyro * d_chi; %7 tesla times 42.6 MHz * 0.1 ppm;
Hct = 0.40; %Hemocrit for inside vessel
BOLD_freq = zeros(1, num_spin);

  function BOLD_outside = bold_outside(Constant, Y, d_vessel, theta, spin_pos, space_center)
   %get the Bold magnetic field inside the vessel

   element1 = sin(theta) ^ 2;
   distance = sqrt((spin_pos(1) - space_center) ^ 2 + (spin_pos(2) - space_center) ^ 2);
   element2 = 2 .* spin_pos(1) / distance - 1;

   BOLD_outside = 0.50 * Constant .* (1 - Y) .* element1 .* element2 .* (0.50 * d_vessel / distance)^2;
  end
  
  function BOLD_inside = bold_inside(Constant, Y, d_vessel, theta, spin_pos, space_center)
    %get the Bold magnetic field outside the vessel
    element1 = 3 .* (cos(theta))^2 - 1;

    BOLD_inside = Hct * Constant .* (1 - Y) .* element1 ./ 6;
    
  end

  for ind_spin = 1: num_spin
    distance = sqrt( (spin_pos(1, ind_spin) - space_center) ^ 2 + (spin_pos(2, ind_spin) - space_center) ^ 2);
      if distance > d_vessel / 2
        BOLD_freq(ind_spin) = bold_outside(Constant, Y, d_vessel, theta, spin_pos(: , ind_spin), space_center);
      else
        BOLD_freq(ind_spin) = bold_inside(Constant, Y, d_vessel, theta, spin_pos(: , ind_spin), space_center);
      end
  end
      
end