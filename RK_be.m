%general RK BE solver
function RKBE_dM = RK_be(mx, my, mz, dwx, dwy, df, T1, T2, dt)
  %Bloch equation of x
  function BE_x = bex(mx, my, mz, dwy, df, T2)
    E1 = my .* df;
    E2 = mz .* dwy;
    E3 = mx ./ T2;
    BE_x = E1 - E2 - E3;
  end
  
  %Bloch Equation of y
  function BE_y = bey(mx, my, mz, dwx, df, T2);
    E1 = mx .* df;
    E2 = mz .* dwx;
    E3 = my ./ T2;
    BE_y = -E1 + E2 - E3;
  end
  
  function BE_z = bez(mx, my, mz, dwx, dwy, T1);
    E1 = mx .* dwy;
    E2 = my .* dwx;
    E3 = (mz - 1) ./ T1;
    BE_z = E1 - E2 - E3;
  end
  
  kx_1 = bex(mx, my, mz, dwy, df, T2);
  ky_1 = bey(mx, my, mz, dwx, df, T2);
  kz_1 = bez(mx, my, mz, dwx, dwy, T1);

  mx_temp = mx + 0.5 .* kx_1 .* dt;
  my_temp = my + 0.5 .* ky_1 .* dt;
  mz_temp = mz + 0.5 .* kz_1 .* dt;

  kx_2 = bex(mx_temp, my_temp, mz_temp, dwy, df, T2);
  ky_2 = bey(mx_temp, my_temp, mz_temp, dwx, df, T2);
  kz_2 = bez(mx_temp, my_temp, mz_temp, dwx, dwy, T1);

  mx_temp = mx + 0.5 .* kx_2 .* dt;
  my_temp = my + 0.5 .* ky_2 .* dt;
  mz_temp = mz + 0.5 .* kz_2 .* dt;

  kx_3 = bex(mx_temp, my_temp, mz_temp, dwy, df, T2);
  ky_3 = bey(mx_temp, my_temp, mz_temp, dwx, df, T2);
  kz_3 = bez(mx_temp, my_temp, mz_temp, dwx, dwy, T1);

  mx_temp = mx + kx_3 .* dt;
  my_temp = my + ky_3 .* dt;
  mz_temp = mz + kz_3 .* dt;

  kx_4 = bex(mx_temp, my_temp, mz_temp, dwy, df, T2);
  ky_4 = bey(mx_temp, my_temp, mz_temp, dwx, df, T2);
  kz_4 = bez(mx_temp, my_temp, mz_temp, dwx, dwy, T1);

  RK_kx = (kx_1 + 2 .* kx_2 + 2 .* kx_3 + kx_4) ./ 6;
  RK_ky = (ky_1 + 2 .* ky_2 + 2 .* ky_3 + ky_4) ./ 6;
  RK_kz = (kz_1 + 2 .* kz_2 + 2 .* kz_3 + kz_4) ./ 6;
  
  RKBE_dM = [RK_kx; RK_ky; RK_kz] .* dt;
end
    