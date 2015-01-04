%Simulate Brain Tumor and healthy tissue under active feedback MRI \
%Using BOLD model to simulate the magnetic field in brain. Using Random walk to simulate the diffusion of protons
%Pulse Sequence used is Spin Echo Sequence

%Programmed by Song Yang 09/01/2014


%Space class to represent a cubic voxel with a cylinder vessel inside
classdef Space<handle
  properties
    %Vessel Size
    d_vessel % vassel diameter um
    S_vessel %Area of vassel in um^2
    
    vessel_ratio%size of vassel ratio comparted to total space
    len_space%the length of cubic space
    
    %coordination of vessel center
    vessel_center
    vessel_angle 
    %fieldmap
    df_var

  end
  
  methods
    function obj = Space(diameter, size_ratio, vessel_angle)
      %initiate parameters
      obj.d_vessel = diameter;
      obj.vessel_ratio = size_ratio;
      obj.vessel_angle = vessel_angle;
      
      obj.S_vessel = pi * (obj.d_vessel / 2) ^ 2;
      obj.len_space = sqrt(obj.S_vessel / obj.vessel_ratio);
      obj.vessel_center = 0.5 * obj.len_space;
    end
     
     %change the vessel angle to B0
    function update_VA(obj)
      obj.vessel_angle = pi * rand();
    end
    
    %plot a intersections of voxel space
    function plot_plane(obj)
      %plot out boundary
      rectangle('Position', [0, 0, obj.len_space, obj.len_space], 'Curvature', [0, 0]);
      hold on;
      %plot inner boundary
      rectangle('Position', [obj.vessel_center - obj.d_vessel / 2, obj.vessel_center - obj.d_vessel / 2, obj.d_vessel, obj.d_vessel], 'Curvature', [1, 1]);
    end
    
    %update df with BOLD model
    function update_bold(obj, pos, Y)
      Nsp = size(pos, 2);
      obj.df_var = zeros(1, Nsp);
      wb = bold_freq(Nsp, obj.d_vessel, obj.vessel_angle, pos, Y, obj.vessel_center);
      obj.df_var = obj.df_var + wb;
    end
    
    %update df with background dipolar field frequency
    function update_dpf(obj, pos, dw_rms)
      dipole_pos = [0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75;
                    0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75;
                    0.25, 0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75] .* obj.len_space;
      r_dipole = 1; % in um
      N_dipole = 8;
      
      Nsp = size(pos, 2);
      dpf = dipole_freq(pos, Nsp, dipole_pos, N_dipole, r_dipole, dw_rms);
      obj.df_var = obj.df_var + dpf;
    end
    
    %plot 3D structure of the model, cube with a vessel inside it
    function plot_3d(obj)
      view(3);
      axis vis3d;
      %plot out boundary
      verts = [0 0 0; 0 1 0; 1 1 0; 1 0 0; 0 0 1; 0 1 1; 1 1 1 ;1 0 1] .* obj.len_space;
      facs = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
      patch('Vertices', verts, 'Faces', facs, 'FaceColor', 'g', 'facealpha', 0.1);
      hold on;
      %plot in boundary
      [X, Y, Z] = cylinder(obj.d_vessel / 2);
      X = X + obj.len_space / 2;
      Y = Y + obj.len_space / 2;
      Z = (Z ./ max(max(Z))) .* obj.len_space; 
      surf(X, Y, Z, 'facealpha', 0.5);
    end
    
    %plot dpf particles inside it
    function plot_dpf(obj)
      [X, Y, Z] = sphere;
      x1 = (X * 0.05 * obj.len_space + 0.25 * obj.len_space);
      x2 = (X * 0.05 * obj.len_space + 0.75 * obj.len_space);
      y1 = (Y * 0.05 * obj.len_space + 0.25 * obj.len_space);
      y2 = (Y * 0.05 * obj.len_space + 0.75 * obj.len_space);
      z1 = (Z * 0.05 * obj.len_space + 0.25 * obj.len_space);
      z2 = (Z * 0.05 * obj.len_space + 0.75 * obj.len_space);
      
      surf(x1, y1, z1);
      hold on;
      surf(x1, y2, z1);
      hold on;
      surf(x1, y1, z2);
      hold on;
      surf(x1, y2, z2);
      hold on;
      surf(x2, y2, z2);
      hold on;
      surf(x2, y1, z2);
      hold on;
      surf(x2, y1, z1);
      hold on;
      surf(x2, y2, z1);
    end
    
    %plot 3D structure of voxel space and dipolar particles inside it. 
    function plot_3d_dpf(obj)
      obj.plot_3d;
      hold on;
      obj.plot_dpf();
    end
    
    %plot magnetic field map induced by BOLD signal source
    function plot_bold(obj, plane_z, Y)
      dx = 0.1;
      x = 0: dx: obj.len_space;
      Nx = length(x);
      Len_x = Nx * Nx;
      [X_axis, Y_axis] = meshgrid(x);
      Z = plane_z * obj.len_space * ones(1, Len_x);
      temp_x = reshape(X_axis, 1, Len_x);
      temp_y = reshape(Y_axis, 1, Len_x);
      pos = [temp_x; temp_y; Z];
      
      obj.update_bold(pos, Y);
      temp_df = reshape(obj.df_var, Nx, Nx);
      
      figure;
      obj.plot_plane();
      hold on;
      C = contourf(X_axis, Y_axis, temp_df);
      clabel(C);
      hold off;
    end
    
    %plot magnetic field map induced by BOLD and dipolar magnetic field 
    function plot_df(obj, plane_z, Y);
      dx = 0.1;
      x = 0: dx: obj.len_space;
      Nx = length(x);
      Len_x = Nx * Nx;
      [X_axis, Y_axis] = meshgrid(x);
      Z = plane_z * obj.len_space * ones(1, Len_x);
      temp_x = reshape(X_axis, 1, Len_x);
      temp_y = reshape(Y_axis, 1, Len_x);
      pos = [temp_x; temp_y; Z];
      
      dw_rms = 2.36;
      %update both Bold field and dpf
      obj.update_bold(pos, Y);
      obj.update_dpf(pos, dw_rms);
      
      temp_df = reshape(obj.df_var, Nx, Nx);
      
      figure;
      obj.plot_plane();
      hold on;
      C = contourf(X_axis, Y_axis, temp_df);
      clabel(C);
      hold off;
    end
    
    %rotate vessel angle by 10 degree.
    function flip_angle(obj)
      obj.vessel_angle = obj.vessel_angle + pi / 18;
    end
    
  end
  
end

