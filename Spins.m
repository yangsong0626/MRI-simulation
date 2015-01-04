%Simulate Brain Tumor and healthy tissue under active feedback MRI \
%Using BOLD model to simulate the magnetic field in brain. Using Random walk to simulate the diffusion of protons
%Pulse Sequence used is Spin Echo Sequence

%Programmed by Song Yang 09/01/2014

classdef Spins<handle
  properties
    num_spin %the number of spins inside the space
    positions %positions of spins
    
    %magnetization of spins
    M
    Mean_Mx
    Mean_My
    Mean_Mz
    Mxy
    
    %Spin physical constants
    T1 %in s 
    T2 % in s
    
    df0
    v_shift  %Frequency shift in Hz 
  end%properties end
  
  methods
    function obj = Spins(N_spins, Freq_shift, T1, T2)
      obj.num_spin = N_spins;
      obj.v_shift = Freq_shift;
      obj.df0 = Freq_shift;
      
      obj.T1 = T1;
      obj.T2 = T2;
    end
    
    %reset accumulated phase shift back to 0;
    
    function init_magnet(obj)
      %initiate the magnetization of spins
      Mx_0 = zeros(1, obj.num_spin);
      My_0 = zeros(1, obj.num_spin);
      Mz_0 = ones(1, obj.num_spin);
      
      obj.M = [Mx_0; My_0; Mz_0];
    end
    
    function update_magnet(obj, dM)
      %update spin magnetization based on increment given
      obj.M = obj.M + dM;
    end
    
    %Calculate the average Mx, My, Mxy, Mz etc
    function update_MeanM(obj)
      obj.Mean_Mx = mean(obj.M(1, :));
      obj.Mean_My = mean(obj.M(2, :));
      obj.Mean_Mz = mean(obj.M(3, :));
      obj.Mxy = sqrt(obj.Mean_Mx ^ 2 + obj.Mean_My ^ 2);
    end
    
    %update the v_shift to giving df
    function update_df(obj, df)
      obj.v_shift = df;
    end
    
    %reset v_shift back to df0;
    function redeem_df(obj)
      obj.v_shift = obj.df0;
    end
  %-----------------------------------------------------------------------------------------------------------  
    %Pulse Sequence
    
    %Spin Echo
    function SE(obj, pulse_order)
      %update spin magnetization to Spin Echo pulse sequence
      if pulse_order == 1
        obj.M = yrot(pi / 2) * obj.M;
        elseif pulse_order == 2
          obj.M = xrot(pi) * obj.M;
        else
          fprintf('There are only 2 pulses in Spin Echo pulse sequence, please enter 1 or 2 for pulse order');
      end
    end  
    
    %PNNP pulse sequence
    function PNNP(obj, pulse_order, pulse_angle)
      %update PNNP pulse sequence
      if pulse_order == 1
        obj.M = xrot(pulse_angle) * obj.M;
        elseif pulse_order == 2
          obj.M = xrot(-pulse_angle) * obj.M;
        elseif pulse_order == 3
          obj.M = yrot(pulse_angle) * obj.M;
        elseif pulse_order == 4
          obj.M = yrot(-pulse_angle) * obj.M;
        else
          fprintf('There are only 4 pulses in PNNP pulse sequence')
      end
    end
    
    %Inversion Recovery
    function IR(obj, pulse_order)
      if pulse_order == 1
        obj.M = xrot(pi) * obj.M;
        elseif pulse_order == 2
          obj.M = xrot(pi / 2) * obj.M;
        elseif pulse_order == 3
          obj.M = yrot(pi) * obj.M;
      end
    end
    
    %Steady State pulse sequence
    function SS(obj, pulse_angle)
      obj.M = xrot(pulse_angle) * obj.M;
    end
    
  %-------------------------------------------------------------------------------------------------------------

  %---------------------------------------------------------------------------------------------------------
  %magnetization evolve
    
    %M evolve under radiation damping field
    function RD_evolve(obj, mxm, mym, tr, fai, dt)
      obj.dwx = (mxm * sin(fai) - mym * cos(fai)) / tr;
      obj.dwy = (mxm * cos(fai) + mym * sin(fai)) / tr;
      dM = RK_be(obj.M(1, :), obj.M(2, :), obj.M(3, :), obj.dwx, obj.dwy, obj.v_shift, obj.T1, obj.T2, dt);
      obj.update_magnet(dM);
    end
    
    %M evolve without radiation damping field
    function FP_evolve(obj, dt)
      dwx = zeros(1, obj.num_spin);
      dwy = zeros(1, obj.num_spin);
      dM = RK_be(obj.M(1, :), obj.M(2, :), obj.M(3, :), dwx, dwy, obj.v_shift, obj.T1, obj.T2, dt);
      obj.update_magnet(dM);
      obj.update_MeanM();
    end
    
    %M evolve without radiation damping field, calculate through rotation Matrix method.
    function FP_RMevolve(obj, dt)
      for isp = 1: obj.num_spin
        rot_mat = freeprecess(dt, obj.T1, obj.T2, obj.v_shift(isp));
        obj.M(:, isp) = rot_mat * obj.M(:, isp);
      end
      obj.update_MeanM();
    end
    
    %M evolving under CW MRI without active feedback field;
    function CW_evolve(obj, cw_freq, dt)
      %obj.M = xrot(cw_angle) * obj.M;
      obj.dwx = cw_freq;
      %obj.dwy = zeros(1, obj.num_spin);
      dM = RK_be(obj.M(1, :), obj.M(2, :), obj.M(3, :), obj.dwx, obj.dwy, obj.v_shift, obj.T1, obj.T2, dt);
      obj.update_magnet(dM);
    end
    
    %M evolve under CW MRI with active feedback field;
    function CWRD_evolve(obj, cw_freq, cw_phi, mxm, mym, tr, rd_fai, dt)
      %obj.M = xrot(cw_angle) * obj.M;
      obj.dwx = (mxm * sin(rd_fai) - mym * cos(rd_fai)) / tr + cw_freq * cos(cw_phi);
      obj.dwy = (mxm * cos(rd_fai) + mym * sin(rd_fai)) / tr + cw_freq * sin(cw_phi);
      dM = RK_be(obj.M(1, :), obj.M(2, :), obj.M(3, :), obj.dwx, obj.dwy, obj.v_shift, obj.T1, obj.T2, dt);
      obj.update_magnet(dM);
    end
    
    %M evolve under CW MRI with active feedback field, calculating through rotation matrix method;
    function CWRD_RMevolve(obj, cw_freq, cw_phi, mxm, mym, tr, rd_fai, dt)
      
      for isp = 1: obj.num_spin
        rot_mat = rdrot(mxm, mym, rd_fai, cw_freq, cw_phi, obj.v_shift(isp), tr, dt);
        obj.M(:, isp) = rot_mat * obj.M(:, isp);
      end
      obj.update_MeanM();
    end
  %------------------------------------------------------------------------------
  %Process code, combination of M evolve and Pulse sequence
  function CWRD_RMprocess_spin(obj, cw_freq, cw_phi, mxm, mym, tr, rd_fai, t_tot, dt)
    evolving_t = 0: dt: t_tot;
    Nt = length(evolving_t);
    
    obj.init_magnet();
    
    
    for it = 1: Nt
      obj.CWRD_RMevolve(obj, cw_freq, cw_phi, mxm, mym, tr, rd_fai, dt);
    end
    
  end
  
  %------------------------------------------------------------------------------
  %plot methods
  %plot spin as a sphere centered at its position with a certain radius
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
    %plot random spin magnetization evolviing
    
  end%methods end
  
end%class end
      