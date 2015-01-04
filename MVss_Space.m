%Simulate Brain Tumor and healthy tissue under active feedback MRI \
%Using BOLD model to simulate the magnetic field in brain. Using Random walk to simulate the diffusion of protons
%Pulse Sequence used is Spin Echo Sequence

%Programmed by Song Yang 09/01/2014


%Space class to represent a cubic voxel with a cylinder vessel inside
classdef MVss_Space<handle
  properties
    N_vessel
    %Vessel Size
    S_vessel %Area of vassel in um^2
    
    vessel_ratio%size of vassel ratio comparted to total space
    len_space%the length of cubic space
    
    %coordination of vessel center
    vessel_strts
    vessel_ends
    vessel_angles 
    vessel_vects
    d_vessels % vassel diameter um
    mean_d
    %fieldmap
    df_var

  end
  
  methods
    function obj = MVss_Space(N_vss, diameter, f_blood, vessel_angles)
      obj.N_vessel = N_vss;
      obj.d_vessels = diameter * ones(1, obj.N_vessel);
      obj.mean_d = mean(obj.d_vessels);
      obj.vessel_ratio = f_blood;
      obj.S_vessel = obj.N_vessel * pi * (obj.mean_d / 2) ^ 2;
      
      obj.len_space = sqrt(obj.S_vessel / obj.vessel_ratio);
      
      obj.vessel_strts = rand(3, N_vss) .* obj.len_space;
      obj.vessel_angles = vessel_angles;
      temp_z = obj.vessel_strts(3, :) + cos(obj.vessel_angles);
      temp_angle = rand(1, obj.N_vessel) .* 2 .* pi;
      temp_x = sin(temp_angle) .* sin(obj.vessel_angles) + obj.vessel_strts(1, :);
      temp_y = cos(temp_angle) .* sin(obj.vessel_angles) + obj.vessel_strts(2, :);
      clear temp_angle;
      obj.vessel_ends = [temp_x; temp_y; temp_z];
      obj.vessel_vects = obj.vessel_strts - obj.vessel_ends;
      %obj.vessel_angles = acos(dot(obj.vessel_vects, z_vect) / (norm(obj.vessel_vects) * norm(z_vect)));
    end
    
    function plot_plane(obj, z_coord)
      %plot out boundary
      rectangle('Position', [0, 0, obj.len_space, obj.len_space], 'Curvature', [0, 0]);
      hold on;
      %plot inner boundary
      %rectangle('Position', [obj.vessel_center - obj.d_vessel / 2, obj.vessel_center - obj.d_vessel / 2, obj.d_vessel, obj.d_vessel], 'Curvature', [1, 1]);
    end
    
    %update df with BOLD model
    function update_bold(obj, pos, dx_B)
      Nsp = size(pos, 2);
      obj.df_var = zeros(1, Nsp);
      for isp = 1: Nsp
        wb = obj.BOLD_freq(pos(:, isp), dx_B);
        obj.df_var(isp) = wb;
      end
    end
    
    function bold_freq = BOLD_freq(obj, pos, dx_B)
      temp_freq = 0;
      gyro = 42.6e6;
      for ivss = 1: obj.N_vessel
        [d_temp, in_temp] = obj.DistanceToVss(ivss, pos);
        if in_temp
          E2 = 3 * (cos(obj.vessel_angles(ivss))) ^ 2 - 1;
          temp_freq = temp_freq + gyro * dx_B * E2;
        else
          if obj.vessel_angles(ivss) ~= 0
            b0_projection = [0; 0; 1] - (dot([0; 0; 1], obj.vessel_vects(:, ivss))) * obj.vessel_vects(:, ivss) ./ norm(obj.vessel_vects(:, ivss));
            temp_vect = pos - obj.vessel_strts(:, ivss);
            pos_projection = (temp_vect) - (dot(temp_vect, obj.vessel_vects(:, ivss))) * obj.vessel_vects(:, ivss) ./ norm(obj.vessel_vects(:, ivss));
            E1 = 2 * (dot(b0_projection, pos_projection)/ (norm(b0_projection) * norm(pos_projection))) ^ 2 - 1;
            E2 = (sin(obj.vessel_angles(ivss))) ^ 2;
            temp_freq = temp_freq + gyro * dx_B * E1 * E2 * (0.5 * obj.d_vessels(ivss) / d_temp) ^ 2;
          end
        end
      end
      bold_freq = temp_freq;
    end

    
    %compute the distance from spin to all vessels
    function [D_vss, IN_vss] = DistanceToVss(obj, ivss, pos)
      D_vss = 2 * norm(cross(obj.vessel_vects(:, ivss), pos - obj.vessel_strts(:, ivss))) / norm(obj.vessel_vects);
      radius = 0.5 * obj.d_vessels(ivss);
      if D_vss < radius
        IN_vss = 1;
      else
        IN_vss = 0;
      end
    end
    
    %test DistanceToVss method
    function result = test_DistanceToVss(obj, ivss)
      [d, in] = obj.DistanceToVss(ivss, obj.vessel_strts(:, ivss));
      if d == 0 && in == 1
        result = true;
      else
        result = false;
      end
    end
    
    
    %update df with background dipolar field frequency
    function update_dpf(obj, pos)
      dipole_pos = [0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75;
                    0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75;
                    0.25, 0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75] .* obj.len_space;
      r_dipole = 1; % in um
      N_dipole = 8;
      
      Nsp = size(pos, 2);
      dpf = dipole_freq(pos, Nsp, dipole_pos, N_dipole, r_dipole);
      obj.df_var = obj.df_var + dpf;
    end
    
    function h = plot_3d(obj)
      h = figure;
      view(3);
      xlim([0, obj.len_space]);
      ylim([0, obj.len_space]);
      zlim([0, obj.len_space]);
      axis vis3d;
      %plot out boundary
      verts = [0 0 0; 0 1 0; 1 1 0; 1 0 0; 0 0 1; 0 1 1; 1 1 1 ;1 0 1] .* obj.len_space;
      facs = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
      patch('Vertices', verts, 'Faces', facs, 'FaceColor', 'g', 'facealpha', 0.1);
      hold on;
      %plot in boundary
      for ivss = 1: obj.N_vessel
        [X, Y, Z] = cylinder(obj.d_vessels(ivss) / 2);
        X = X + obj.vessel_strts(1, ivss);
        Y = Y + obj.vessel_strts(2, ivss);
        Z = (Z ./ max(max(Z))) .* obj.len_space; 
        vss_plot = surf(X, Y, Z, 'facealpha', 0.5);
        if obj.vessel_angles(ivss) ~= 0
          axis_rot = cross([0; 0; 1], obj.vessel_vects(:, ivss));
          %axis_rot = cross([0; 0; 1], obj.vessel_vects(:, ivss));
          rotate(vss_plot, axis_rot, obj.vessel_angles(ivss) * 180 / pi, obj.vessel_strts(:, ivss));
        end
        hold on;
      end
      hold off;
    end
    
    function h = plot_df(obj, plane_z, Y)
      dx = 1;
      x = 0: dx : obj.len_space;
      Nx = length(x);
      [X_axis, Y_axis] = meshgrid(x);
      Z = plane_z * obj.len_space * ones(Nx, Nx);
      bold_df = zeros(Nx, Nx);
      
      for ix = 1 : Nx
        for iy = 1: Nx
          temp_pos = [X_axis(ix, iy); Y_axis(ix, iy); Z(ix, iy)];
          bold_df(ix, iy) = obj.BOLD_freq(temp_pos, Y);
        end
      end
      h = figure;
      C = contourf(X_axis, Y_axis, bold_df);
      clabel(C);
    end
    
    
  end
  
end

