%===============================================================================================
%Simulate a voxel of Brain Tumor or healthy tissue under MRI 
%Distribute a certain number of spins in a cubic space. @init_magnet() initiate the spins magnetization to (0, 0, 1).
%@init_positions() method assign random positions to spins as their starting points.
%Spins diffusion are models as random walk via method @RW(), and plotted by @plot_RW() method
%Spins evolution are completed by quoting methods from Spins superclass.
%Magnetic field inhomogeneity are calculated by methods from Space superclass. 
%Programmed by Song Yang 09/01/2014
%===============================================================================================================

classdef VXL<Space & Spins
  properties
    %diffusion Constants
    D %Diffusion constant, in um^2/ms
    q = 20; %in /ms
    step_length %random walk step length
    
  end%properties end
  
  methods
    function obj = VXL(N_spins, Freq_shift, T1, T2, diameter, size_ratio, Diffusion, theta)
      %Initiate Spins and Space superclass parameters
      obj = obj@Spins(N_spins, Freq_shift, T1, T2);
      obj = obj@Space(diameter, size_ratio, theta);
      
      %Initiate VXL parameters
      obj.D = Diffusion;
      obj.step_length = sqrt(6 * obj.D / obj.q);
    end
    
    %initiate the positions of spins, within the obj.len_space, and inside or outside the inner boundary
    function init_positions(obj)
      spin_positions = zeros(3, obj.num_spin);
      for idx_spin = 1 : obj.num_spin
        temp_spin_pos = obj.len_space * rand(3, 1);
        temp_distance = sqrt((temp_spin_pos(1) - obj.vessel_center) ^ 2 + (temp_spin_pos(2) - obj.vessel_center) ^ 2);
        
        while temp_distance == obj.d_vessel
          temp_spin_pos = obj.len_space * rand(3, 1);
          temp_distance = sqrt((temp_spin_pos(1) - obj.vessel_center) ^ 2 + (temp_spin_pos(2) - obj.vessel_center) ^ 2);
        end
        
        spin_positions(: , idx_spin) = temp_spin_pos;
      end%for loop end
      
      obj.positions = spin_positions;
    end
    
    %simulate diffusion of spins inside vessel centre using random walk model
    function RW(obj, dt)
      %for spins inside obj.d_vessel (diameter of the cylinder), stay inside. For outside spins, stay outside
      obj.step_length = sqrt(2 * obj.D * dt);
      for idx_spin = 1 : obj.num_spin
        angle_1 = rand * pi;
        angle_2 = rand * 2.0 * pi;
        step_matrix = [sin(angle_1) * cos(angle_2); sin(angle_1) * sin(angle_2); cos(angle_1)] .* obj.step_length;
        
        origin_position = obj.positions(:, idx_spin);
        origin_distance = sqrt((obj.positions(1, idx_spin) - obj.vessel_center) ^ 2 + (obj.positions(2, idx_spin) - obj.vessel_center) ^ 2);
        obj.positions(:, idx_spin) = obj.positions(:, idx_spin) + step_matrix;
        after_distance = sqrt((obj.positions(1, idx_spin) - obj.vessel_center) ^ 2 + (obj.positions(2, idx_spin) - obj.vessel_center) ^ 2);
        
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
        if (origin_distance > obj.d_vessel / 2 && after_distance <= obj.d_vessel / 2) || (origin_distance < obj.d_vessel / 2 && after_distance >= obj.d_vessel / 2)
          obj.positions(1, idx_spin) = 2 * obj.vessel_center - origin_position(1) + step_matrix(1); %obj.positions(1, idx_spin) - step_matrix(1);
          obj.positions(2, idx_spin) = 2 * obj.vessel_center - origin_position(2) + step_matrix(2);%obj.positions(2, idx_spin) - step_matrix(2);
         
        end%end of inner boundary conditions
      end %end of idx_spin for loop
      
    end 
    
    %Plot the random walk track of two spins, one inside obj.d_vessel and the other outside
    function h = plot_RW(obj, time_step, tot_time)
      spin_1 = 1;
      spin_2 = 2;
      spin_3 = 3;
      spin_4 = 4;
      
      obj.init_positions();
      h = figure;
      obj.plot_plane();
      hold on;
      t = 0 : time_step : tot_time;
      num_step = length(t);
      for idx_step = 1 : num_step
        obj.RW(time_step);
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
    %initiate spins magnetization
    function init_magnet(obj)
      init_magnet@Spins(obj);
    end
    
    %update spins magnetization by adding increment dM 1*3 matrix.
    function update_magnet(obj, dM)
      update_magnet@Spins(obj, dM);
    end
    
    %calculate Mean_Mx, Mean_My, Mean_Mz, Mxy
    function update_MeanM(obj)
      update_MeanM@Spins(obj);
    end
    
    %update spins off-resonance frequency by adding increment df
    function update_df(obj, df)
      update_df@Spins(obj, df);
    end
    
    %set spins off-resonance frequency back to initial value
    function redeem_df(obj)
      redeem_df@Spins(obj);
    end
    
    %update spins off-resonance frequency by adding spins and space inhomogeneity up. 
    function update_spdf(obj)
      obj.v_shift = obj.v_shift + obj.df_var;
    end
  
    %pulse sequences
    %Spin Echo pulse sequence
    function SE(obj, pls_ord)
      SE@Spins(obj, pls_ord);
    end
    
    %Postive negative negative positive pulse sequence
    function PNNP(obj, pls_ord, pls_angle)
      PNNP@Spins(obj, pls_ord, pls_angle);
    end
    
    %inversion reversion pulse sequence
    function IR(obj, pls_ord)
      IR@Spins(obj, pls_or);
    end
  
    %Steady state pulse sequence
    function SS(obj, pls_angle)
      SS@Spins(obj, pls_angle);
    end
    
    
    %M evolve
    %spin evolution under radiation damping field
    function RD_evolve(obj, mxm, mym, tr, fai, dt)
      RD_evolve@Spins(obj, mxm, mym, tr, fai, dt);
    end
    
    %spin evolution by free precessing
    function FP_evolve(obj, dt)
      FP_evolve@Spins(obj, dt);
    end
    
    %spin evolution under continuous wave field
    function CW_evolve(obj, cw_angle, dt)
      CW_evolve@Spins(obj, cw_angle, dt);
    end
    
    %spin evolution under CW field and radiation damping field.
    function CWRD_evolve(obj, cw_angle, mxm, mym, tr, rd_fai, dt)
      CWRD_evolve@Spins(obj, cw_angle, mxm, mym, tr, rd_fai, dt);
    end
    
    %Spin evolution under CW field and RD field via rotating matrix method.
    function CWRD_RMevolve(obj, cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt)
      CWRD_RMevolve@Spins(obj, cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
    end
    
    %single Voxel CWRD_RM process in t_tot time, with time step dt.
    function CWRD_RMprocess(obj, cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y_vx, rotating_angle)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      
      for it = 1: Nt
        obj.RW(dt);
        obj.update_bold(Y_vx);
        obj.update_spdf();
        mxm = obj.Mean_Mx;
        mym = obj.Mean_My;
        obj.CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
        obj.update_MeanM();
        if rotating_angle
          obj.flip_angle();
        end
        
        obj.redeem_df();
      end
    end
    
    %test CWRD process
    function h = test_CWprocess(obj, cw_angle, cw_phi, tr, rd_fai, t_tot, dt, Y_vx, rotating_angle)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      
      Mzt = zeros(1, Nt);
      Mxyt = zeros(1, Nt);
      for it = 1: Nt
        obj.RW(dt);
        obj.update_bold(Y_vx);
        obj.update_spdf();
        mxm = obj.Mean_Mx;
        mym = obj.Mean_My;
        obj.CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
        obj.update_MeanM();
        if rotating_angle
          obj.flip_angle();
        end
        
        obj.redeem_df();
        Mzt(it) = obj.Mean_Mz;
        Mxyt(it) = obj.Mxy;
      end
      
      h = figure;
      plot(evolving_t, Mzt, 'g', evolving_t, Mxyt, 'b');
      legend('Mz', 'Mxy');
      
      title('Signal evolving with time under CW MRI, with BOLD perturbation only');
      xlabel('t(s)');
      ylabel('M');
    end
    
    %Spin Echo evolution process of a total time length of t_tot, with time interval of dt.
    function SE_process(obj, N_te, Y_vx, t_tot, dt, rotating_angle)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      obj.SE(1);
      
      for it = 1: Nt
        obj.RW(dt);
        if it == N_te
          obj.SE(2);
        end
        obj.update_bold(Y_vx);
        obj.update_spdf();
        obj.FP_evolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
        if rotating_angle
          obj.flip_angle();
        end
      end
      
    end
    
    %test SE_process
    function h = test_SEprocess(obj, N_te, Y_vx, t_tot, dt, rotating_angle)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      Mxyt = zeros(1, Nt);
      Mzt = zeros(1, Nt);
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      obj.SE(1);
      
      for it = 1: Nt
        obj.RW(dt);
        if it == N_te
          obj.SE(2);
        end
        obj.update_bold(Y_vx);
        obj.update_spdf();
        obj.FP_evolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
        if rotating_angle
          obj.flip_angle();
        end
        Mzt(it) = obj.Mean_Mz;
        Mxyt(it) = obj.Mxy;
      end
      h = figure;
      plot(evolving_t, Mzt, 'g', evolving_t, Mxyt, 'b');
      legend('Mz', 'Mxy');
      
      title('Signal evolving with time under SE MRI, with BOLD perturbation only');
      xlabel('t(s)');
      ylabel('M');
    end
    
    %Free precessing process
    function FP_process(obj, N_te, Y_vx, t_tot, dt, rotating_angle)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      
      for it = 1: Nt
        obj.RW(dt);
        obj.update_bold(Y_vx);
        obj.update_spdf();
        obj.FP_evolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
        if rotating_angle
          obj.flip_angle();
        end
      end
      
    end
    
    %plot spins in the space, for demonstration only!
    function plot_spin(obj, radius, number)
      [X, Y, Z] = sphere;
      for isp = 1: number
        x_temp = X * radius + obj.positions(1, isp);
        y_temp = Y * radius + obj.positions(2, isp);
        z_temp = Z * radius + obj.positions(3, isp);
        surf(x_temp, y_temp, z_temp);
        hold on;
      end
      hold off;
    end
    %--------------------------------------------------------------------------
    %Space methods
    function plot_plane(obj)
      plot_plane@Space(obj);
    end
    
    function plot_3d(obj)
      plot_3d@Space(obj);
    end
    
    function plot_3d_dpf(obj)
      plot_3d_dpf@Space(obj);
    end
    
    function update_bold(obj, Y)
      update_bold@Space(obj, obj.positions, Y);
    end
    
    function update_dpf(obj)
      update_dpf@Space(obj, obj.positions);
    end
    
    function update_VA(obj)
      update_VA@Space(obj);
    end
    
    %VXL methods
    function flip_angle(obj)
      flip_angle@Space(obj);
    end
  end%methods end
end%class end
      