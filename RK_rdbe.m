function [RK_kx, RK_ky, RK_kz] = RK_rdbe(mx, my, mz, mxm, mym, tr, fai, dw, dt)
%Runge_kutta method to solve Bloch equation of Mx
  function BE_x = bex(my, mz, mxm, tr, dw)
  %bloch equation of MX
  E2 = my .* dw;
  E3 = mz .* mxm ./ tr;

  BE_x = E2 - E3;
  end
  
  function BE_y = bey(mx, mz, mym, tr, dw)
  %Bloch Equation of My
  E2 = mx .* dw;
  E3 = mz .* mym ./ tr;

  BE_y = - E2 - E3;
  end
  
  function BE_z = bez(mx, my, tr, mxm, mym)
  %Bloch equation of Mz
  element1 = (mx .* mxm + my .* mym);
  BE_z = element1 ./ tr;
  end
  
kx_1 = bex(my, mz, mxm, tr, dw);
ky_1 = bey(mx, mz, mym, tr, dw);
kz_1 = bez(mx, my, tr, mxm, mym);

mx_temp = mx + 0.5 .* kx_1 .* dt;
my_temp = my + 0.5 .* ky_1 .* dt;
mz_temp = mz + 0.5 .* kz_1 .* dt;

kx_2 = bex(my_temp, mz_temp, mxm, tr, dw)
ky_2 = bey(mx_temp, mz_temp, mym, tr, dw);
kz_2 = bez(mx_temp, my_temp, tr, mxm, mym);

mx_temp = mx + 0.5 .* kx_2 .* dt;
my_temp = my + 0.5 .* ky_2 .* dt;
mz_temp = mz + 0.5 .* kz_2 .* dt;

kx_3 = bex(my_temp, mz_temp, mxm, tr, dw);
ky_3 = bey(mx_temp, mz_temp, mym, tr, dw);
kz_3 = bez(mx_temp, my_temp, tr, mxm, mym);

mx_temp = mx + kx_3 .* dt;
my_temp = my + ky_3 .* dt;
mz_temp = mz + kz_3 .* dt;

kx_4 = bex(my_temp, mz_temp, mxm, tr, dw);
ky_4 = bey(mx_temp, mz_temp, mym, tr, dw);
kz_4 = bez(mx_temp, my_temp, tr, mxm, mym);

RK_kx = (kx_1 + 2 .* kx_2 + 2 .* kx_3 + kx_4) ./ 6;
RK_ky = (ky_1 + 2 .* ky_2 + 2 .* ky_3 + ky_4) ./ 6;
RK_kz = (kz_1 + 2 .* kz_2 + 2 .* kz_3 + kz_4) ./ 6;
end
  