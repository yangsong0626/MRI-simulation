%compute the rotation matrix for RD field
function RDrot = rdrot(mxm, mym, fai, df, df_phi, dw, tr, dt)
  x1 = df * cos(df_phi) + (mxm * sin(fai) - mym * cos(fai)) / (2 * pi * tr);
  y1 = df * sin(df_phi) + (mxm * cos(fai) + mym * sin(fai)) / (2 * pi * tr);
  z1 = dw;
  
  B1 = [x1; y1; z1];
  B1_strength = norm(B1);
  
 angle1 = -2 * pi * B1_strength * dt;
  x1 = x1 / B1_strength;
  y1 = y1 / B1_strength;
  z1 = z1 / B1_strength;
  
  RDrot = [cos(angle1)+(1-cos(angle1))*x1^2 (1-cos(angle1))*x1*y1-z1*sin(angle1) (1-cos(angle1))*x1*z1+y1*sin(angle1);(1-cos(angle1))*y1*x1+z1*sin(angle1) cos(angle1)+(1-cos(angle1))*y1^2 (1-cos(angle1))*y1*z1-x1*sin(angle1);(1-cos(angle1))*z1*x1-y1*sin(angle1) (1-cos(angle1))*z1*y1+x1*sin(angle1) cos(angle1)+(1-cos(angle1))*z1^2];