%Solve Bloch Equation of Free Precess using Ruger Kutta Methods
function RK_dM = RK_fpbe(mx, my, mz, df, dt, T1, T2)
  %Bloch Equation of Mx
  function BE_x = bex(mx, my, df, T2)
    E1 = df .* my;
    E2 = mx ./ T2;
    BE_x = E1 - E2;
  end
  
  %Bloch Equation of My
  function BE_y = bey(mx, my, df, T2)
    E1 = -df .* mx;
    E2 = -my ./ T2;
    BE_y = E1 + E2;
  end
  
  %Bloch Euqation of Mz
  function BE_z = bez(mz, T1)
    BE_z = -(mz - 1) ./ T1;
  end
  
  kx_1 = bex(mx, my, df, T2);
  ky_1 = bey(mx, my, df, T2);
  kz_1 = bez(mz, T1);
  
  mx_temp = mx + 0.5 .* kx_1 .* dt;
  my_temp = my + 0.5 .* ky_1 .* dt;
  mz_temp = mz + 0.5 .* kz_1 .* dt;

  kx_2 = bex(mx_temp, my_temp, df, T2);
  ky_2 = bey(mx_temp, my_temp, df, T2);
  kz_2 = bez(mz_temp, T1);

  mx_temp = mx + 0.5 .* kx_2 .* dt;
  my_temp = my + 0.5 .* ky_2 .* dt;
  mz_temp = mz + 0.5 .* kz_2 .* dt;

  kx_3 = bex(mx_temp, my_temp, df, T2);
  ky_3 = bey(mx_temp, my_temp, df, T2);
  kz_3 =  bez(mz_temp, T1);

  mx_temp = mx + kx_3 .* dt;
  my_temp = my + ky_3 .* dt;
  mz_temp = mz + kz_3 .* dt;

  kx_4 = bex(mx_temp, my_temp, df, T2);
  ky_4 = bey(mx_temp, my_temp, df, T2);
  kz_4 = bez(mz_temp, T1);

  RK_kx = (kx_1 + 2 .* kx_2 + 2 .* kx_3 + kx_4) ./ 6;
  RK_ky = (ky_1 + 2 .* ky_2 + 2 .* ky_3 + ky_4) ./ 6;
  RK_kz = (kz_1 + 2 .* kz_2 + 2 .* kz_3 + kz_4) ./ 6;

  RK_dM = [RK_kx; RK_ky; RK_kz] .* dt;
end