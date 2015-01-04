%Simulate Brain Tumor and healthy tissue under active feedback MRI 
%Using BOLD model to simulate the magnetic field in brain. Using Random walk to simulate the diffusion of protons
%Pulse Sequence used is Spin Echo Sequence

%Programmed by Song Yang 09/01/2014

classdef Mvv_VXL<MVss_Space & Spins
  properties
    %diffusion Constants
    D %Diffusion constant, in um^2/ms
    q = 20; %in /ms
    step_length %random walk step length
  
    
  end%properties end
  
  methods
    function obj = Mvv_VXL(N_spins, Freq_shift, T1, T2, N_vss, d_vss, f_blood, vessel_angles, Diffusion)
      %Vxl
      obj = obj@Spins(N_spins, Freq_shift, T1, T2);
      obj = obj@MVss_Space(N_vss, d_vss, f_blood, vessel_angles);
      
      obj.D = Diffusion;
    end
    
    %initiate the positions of spins, within the obj.len_space, and inside or outside the inner boundary
    function init_positions(obj)
      spin_positions = zeros(3, obj.num_spin);
      for idx_spin = 1 : obj.num_spin
        temp_spin_pos = obj.len_space * rand(3, 1);
        for ivss = 1: obj.N_vessel
          [temp_distance, temp_in] = obj.DistanceToVss(ivss, temp_spin_pos);
          while temp_distance == obj.d_vessels(ivss)
            temp_spin_pos = obj.len_space * rand(3, 1);
            [temp_distance, temp_in] = obj.DistanceToVss(ivss, temp_spin_pos);
          end
        end
        spin_positions(: , idx_spin) = temp_spin_pos;
      end%for loop end
      
      obj.positions = spin_positions;
    end
    
    %simulate diffusion of spins inside vessel center using random walk model
    function RW(obj, dt)
      %for spins inside obj.d_vessel (diameter of the cylinder), stay inside. For outside spins, stay outside
      obj.step_length = sqrt(2 * obj.D * dt);
      for idx_spin = 1 : obj.num_spin
        angle_1 = rand * pi;
        angle_2 = rand * 2.0 * pi;
        step_matrix = [sin(angle_1) * cos(angle_2); sin(angle_1) * sin(angle_2); cos(angle_1)] .* obj.step_length;
        
        origin_position = obj.positions(:, idx_spin);
        obj.positions(:, idx_spin) = obj.positions(:, idx_spin) + step_matrix;
        
        %Out boundary conditions
        for dummy_idx = 1 : 3
          if obj.positions(dummy_idx, idx_spin) > obj.len_space
            obj.positions(dummy_idx, idx_spin) = obj.positions(dummy_idx, idx_spin) - obj.len_space;
          end
          
          if obj.positions(dummy_idx, idx_spin) < 0
            obj.positions(dummy_idx, idx_spin) = obj.positions(dummy_idx, idx_spin) + obj.len_space;
          end
        end%for loop end
        
        %Inner boundary conditions
        for ivss = 1 : obj.N_vessel
          [origin_distance, origin_in] = obj.DistanceToVss(ivss, origin_position);
          [after_distance, after_in] = obj.DistanceToVss(ivss, obj.positions(:, idx_spin));
          if (origin_in == 1 && after_in == 0) || (origin_in == 0 && after_in == 1)
            obj.positions(:, idx_spin) = origin_position;
          end
        
        end%end of inner boundary conditions
      end %end of idx_spin for loop
      
    end 
    
    %Plot the random walk track of two spins, one inside obj.d_vessel and the other outside
    function plot_RW(obj, time_step, tot_time)
      spin_1 = 1;
      spin_2 = 2;
      spin_3 = 3;
      spin_4 = 4;
      
      figure(1);
      hold on;
      t = 0 : time_step : tot_time;
      num_step = length(t);
      for idx_step = 1 : num_step
        obj.RW();
        scatter(obj.positions(1, spin_1), obj.positions(2, spin_1), 20, 'r', 'o');
        scatter(obj.positions(1, spin_2), obj.positions(2, spin_2), 20, 'b', '+');
        scatter(obj.positions(1, spin_3), obj.positions(2, spin_3), 20, 'y', '*');
        scatter(obj.positions(1, spin_4), obj.positions(2, spin_4), 20, 'g', 'X');
        legend('Spin 1', 'Spin 2', 'Spin 3', 'Spin 4');
        hold on;
      end
      hold off;
    end
    
    %------------------------------------------------------------------------------
    %spins methods
    function init_magnet(obj)
      init_magnet@Spins(obj);
    end
    
    function update_magnet(obj, dM)
      update_magnet@Spins(obj, dM);
    end
    
    function update_MeanM(obj)
      update_MeanM@Spins(obj);
    end
    
    function update_df(obj, df)
      update_df@Spins(obj, df);
    end
    
    function redeem_df(obj)
      redeem_df@Spins(obj);
    end
    
    function update_spdf(obj)
      obj.v_shift = obj.v_shift + obj.df_var;
    end
  
    %pulse sequences
    function SE(obj, pls_ord)
      SE@Spins(obj, pls_ord);
    end
    
    function PNNP(obj, pls_ord, pls_angle)
      PNNP@Spins(obj, pls_ord, pls_angle);
    end
    
    function IR(obj, pls_ord)
      IR@Spins(obj, pls_or);
    end
  
    
    function SS(obj, pls_angle)
      SS@Spins(obj, pls_angle);
    end
    
    
    %M evolve
    function RD_evolve(obj, mxm, mym, tr, fai, dt)
      RD_evolve@Spins(obj, mxm, mym, tr, fai, dt);
    end
    
    function FP_evolve(obj, dt)
      FP_evolve@Spins(obj, dt);
    end
    
    function FP_RMevolve(obj, dt)
      FP_RMevolve@Spins(obj, dt);
    end
    
    function CW_evolve(obj, cw_angle, dt)
      CW_evolve@Spins(obj, cw_angle, dt);
    end
    
    function CWRD_evolve(obj, cw_angle, mxm, mym, tr, rd_fai, dt)
      CWRD_evolve@Spins(obj, cw_angle, mxm, mym, tr, rd_fai, dt);
    end
    
    function CWRD_RMevolve(obj, cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt)
      CWRD_RMevolve@Spins(obj, cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
    end
    
    %--------------------------------------------------------------------------
    %Space methods
    function plot_plane(obj)
      plot_plane@MVss_Space(obj);
    end
    
    function plot_3d(obj)
      plot_3d@MVss_Space(obj);
    end
    
    function update_bold(obj, dxB)
      update_bold@MVss_Space(obj, obj.positions, dxB);
    end
    
    function update_dpf(obj)
      update_dpf@MVss_Space(obj, obj.positions);
    end
    
    function update_VA(obj)
      update_VA@MVss_Space(obj);
    end
    
    function [D_vss, IN_vss] = DistanceToVss(obj, ivss, pos)
      [D_vss, IN_vss] = DistanceToVss@MVss_Space(obj, ivss, pos);
    end
    
    function h = plot_df(obj, plane_z, Y)
      h = plot_df@MVss_Space(obj, plane_z, Y);
    end
    
    %VXL methods
    function flip_angle(obj)
      flip_angle@Space(obj);
    end
    
    function FP_process(obj, Y_saturation, dt, t_tot)
      Nt = t_tot / dt + 1;
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      for it = 1:Nt
        obj.RW(dt);
        obj.update_bold(Y_saturation);
        obj.update_spdf();
        obj.FP_RMevolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
      end
    end
    
    function CWRD_process(obj, Y_saturation, cw_angle, cw_phi, tr, phas, dt, t_tot)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      for it = 1: Nt
        mxm = obj.Mean_Mx;
        mym = obj.Mean_My;
        obj.RW(dt);
        obj.update_bold(Y_saturation);
        obj.update_spdf();
        obj.CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, phas, dt);
        obj.update_MeanM();
        obj.redeem_df();
      end
    end
    function h = test_CWRDprocess(obj, Y_saturation, cw_angle, cw_phi, tr, phas, dt, t_tot)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      Mxyt = zeros(1, Nt);
      Mzt = zeros(1, Nt);
      
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      for it = 1: Nt
        mxm = obj.Mean_Mx;
        mym = obj.Mean_My;
        obj.RW(dt);
        obj.update_bold(Y_saturation);
        obj.update_spdf();
        obj.CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, phas, dt);
        obj.update_MeanM();
        obj.redeem_df();
        Mxyt(it) = obj.Mxy;
        Mzt(it) = obj.Mean_Mz;
      end
      h = figure;
      plot(evolving_t, Mxyt, 'g', evolving_t, Mzt, 'b');
      legend('Mxy', 'Mz');
      xlabel('t(s)');
      ylabel('M');
      title('Test Single Mvv VXL CW process');
    end
      
    function SE_process(obj, Y_saturation, dt, t_tot, Nte)
      Nt = t_tot / dt + 1;
      obj.init_positions();
      obj.init_magnet();
      obj.SE(1);
      obj.update_MeanM();
      for it = 1:Nt
        obj.RW(dt);
        obj.update_bold(Y_saturation);
        obj.update_spdf();
        if it == Nte
          obj.SE(2);
        end
        obj.FP_RMevolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
      end
    end
    
    function h = test_SEprocess(obj, Y_saturation, dt, t_tot, Nte)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      Mxyt = zeros(1, Nt);
      
      obj.init_positions();
      obj.init_magnet();
      obj.SE(1);
      obj.update_MeanM();
      for it = 1:Nt
        obj.RW(dt);
        obj.update_bold(Y_saturation);
        obj.update_spdf();
        if it == Nte
          obj.SE(2);
        end
        obj.FP_RMevolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
        Mxyt(it) = obj.Mxy;
      end
      
      h = figure;
      plot(evolving_t, Mxyt, 'b');
      xlabel('t(s)');
      ylabel('M');
      title('test Single Mvv VXL SE process');
    end
      
  end%methods end
end%class end
      