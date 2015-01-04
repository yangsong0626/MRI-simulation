%Voxel group is a set of voxels with different properties

classdef TwoVoxelGroup<handle
  properties
    %Voxels
    voxels_1
    voxels_2
    r12
    N_voxel
    
    %Voxel posistion

    voxel_angle
    
    %Magnetization signal
    Sx_1
    Sx_2
    Sy_1
    Sy_2
    Sz_1
    Sz_2
    Sxy_1
    Sxy_2
    
    %oxygen saturation level
    Y_vxg1
    Y_vxg2
  end%properties end
  
  methods
    function obj = TwoVoxelGroup(N_voxel, num_spin, v, T1, T2, diameter, rsize, Y, Diffusion, theta, rvxg)
      obj.voxels_1 = VXL.empty(0, N_voxel);
      obj.voxels_2 = VXL.empty(0, N_voxel);
      obj.r12 = rvxg;
      for ivx = 1: N_voxel
          obj.voxels_1(ivx) = VXL(num_spin, v(1), T1, T2, diameter, rsize(1), Diffusion, theta(ivx));
          obj.voxels_2(ivx) = VXL(num_spin, v(2), T1, T2, diameter, rsize(2), Diffusion, theta(ivx));
      end
      obj.N_voxel = N_voxel;
      
      obj.Y_vxg1 = Y(1);
      obj.Y_vxg2 = Y(2);
    end
    
    function VGMean(obj)
      obj.Sx_1 = 0;
      obj.Sx_2 = 0;
      obj.Sy_1 = 0;
      obj.Sy_2 = 0;
      obj.Sz_1 = 0;
      obj.Sz_2 = 0;
      obj.Sxy_1 = 0;
      obj.Sxy_2 = 0;
      for ivx = 1: obj.N_voxel
        obj.Sx_1 = obj.Sx_1 + obj.voxels_1(ivx).Mean_Mx / obj.N_voxel;
        obj.Sx_2 = obj.Sx_2 + obj.voxels_2(ivx).Mean_Mx / obj.N_voxel;
        obj.Sy_1 = obj.Sy_1 + obj.voxels_1(ivx).Mean_My / obj.N_voxel;
        obj.Sy_2 = obj.Sy_2 + obj.voxels_2(ivx).Mean_My / obj.N_voxel;
        obj.Sxy_1 = obj.Sxy_1 + obj.voxels_1(ivx).Mxy / obj.N_voxel;
        obj.Sxy_2 = obj.Sxy_2 + obj.voxels_2(ivx).Mxy / obj.N_voxel;
        obj.Sz_1 = obj.Sz_1 + obj.voxels_1(ivx).Mean_Mz / obj.N_voxel;
        obj.Sz_2 = obj.Sz_2 + obj.voxels_2(ivx).Mean_Mz / obj.N_voxel;
      end
    end
    
    function init_magnet(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).init_magnet();
        obj.voxels_2(ivx).init_magnet();
      end
    end
    
    function init_positions(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).init_positions();
        obj.voxels_2(ivx).init_positions();
      end
    end
    
    function RW(obj, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).RW(dt);
        obj.voxels_2(ivx).RW(dt);
      end
    end
    
    function update_magnet(obj, dM)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).update_magnet(dM);
        obj.voxels_2(ivx).update_magnet(dM);
      end
    end
    
    function update_MeanM(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).update_MeanM();
        obj.voxels_2(ivx).update_MeanM();
      end
    end
    
    function update_df(obj, df)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).update_df(df);
        obj.voxels_2(ivx).update_df(df);
      end
    end
    
    function SE(obj, pls_ord)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).SE(pls_ord);
        obj.voxels_2(ivx).SE(pls_ord);
      end
    end
    
    function PNNP(obj, pls_ord, pls_angle)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).PNNP(pls_ord, pls_angle);
        obj.voxels_2(ivx).PNNP(pls_ord, pls_angle);
      end
    end
    
    function IR(obj, pls_ord)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).IR(pls_ord);
        obj.voxels_2(ivx).IR(pls_ord);
      end
    end
    
    function SS(obj, pls_angle)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).SS(pls_angle);
        obj.voxels_2(ivx).SS(pls_angle);
      end
    end
    
    function RD_evolve(obj, mxm, mym, tr, fai, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).RD_evolve(mxm, mym, tr, fai, dt);
        obj.voxels_2(ivx).RD_evolve(mxm, mym, tr, fai, dt);
      end
    end
    
    function FP_evolve(obj, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).FP_evolve(dt);
        obj.voxels_2(ivx).FP_evolve(dt);
      end
    end
    
    function CW_evolve(obj, cw_angle, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).CW_evolve(cw_angle, dt);
        obj.voxels_2(ivx).CW_evolve(cw_angle, dt);
      end
    end
    
    function CWRD_evolve(obj, cw_angle, mxm, mym, tr, rd_fai, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).CWRD_evolve(cw_angle, mxm, mym, tr, rd_fai, dt);
        obj.voxels_2(ivx).CWRD_evolve(cw_angle, mxm, mym, tr, rd_fai, dt);
      end
    end
    
    function CWRD_RMevolve(obj, cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt)
      for ivx = 1: obj.N_voxel
        obj.voxels_1(ivx).CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
        obj.voxels_2(ivx).CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
      end
    end
    
    function CWRD_RMprocess(obj, cw_angle, cw_phi, tr, rd_fai, t_tot, dt, rotating_angle);
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      obj.VGMean();
      for it = 1: Nt
        mxm_1 = obj.Sx_1; mym_1 = obj.Sy_1;
        mxm_2 = obj.Sx_2; mym_2 = obj.Sy_2;
        mxm = mxm_1 * obj.r12 + mym_2 * (1 - obj.r12);
        mym = mym_1 * obj.r12 + mym_2 * (1 - obj.r12);
        obj.RW(dt);
        obj.update_bold();
        obj.update_spdf();
        
        obj.CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
        obj.update_MeanM();
        obj.redeem_df();
        obj.VGMean();
        obj.flip_angle(rotating_angle);
      end
    end
    
    function h = test_CWprocess(obj, cw_angle, cw_phi, tr, rd_fai, t_tot, dt, rotating_angle);
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      obj.VGMean();
      Mxyt_vxg1 = zeros(1, Nt);
      Mxyt_vxg2 = zeros(1, Nt);
      Mzt_vxg1 = zeros(1, Nt);
      Mzt_vxg2 = zeros(1, Nt);
      for it = 1: Nt
        mxm_1 = obj.Sx_1; mym_1 = obj.Sy_1;
        mxm_2 = obj.Sx_2; mym_2 = obj.Sy_2;
        mxm = mxm_1 * obj.r12 + mym_2 * (1 - obj.r12);
        mym = mym_1 * obj.r12 + mym_2 * (1 - obj.r12);
        obj.RW(dt);
        obj.update_bold();
        obj.update_spdf();
        
        obj.CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
        obj.update_MeanM();
        obj.redeem_df();
        obj.VGMean();
        obj.flip_angle(rotating_angle);
        Mxyt_vxg1(it) = obj.Sxy_1;
        Mxyt_vxg2(it) = obj.Sxy_2;
        Mzt_vxg1(it) = obj.Sz_1;
        Mzt_vxg2(it) = obj.Sz_2;
      end
      contrast_z = abs(Mzt_vxg1 - Mzt_vxg2);
      contrast_xy = abs(Mxyt_vxg1 - Mxyt_vxg2);
      h = figure;
      title('Signal Evolving with time under CW MRI, with Bold only.');
      subplot(121);
      plot(evolving_t, Mxyt_vxg1, 'g', evolving_t, Mxyt_vxg2, 'b', evolving_t, contrast_xy, 'k');
      legend('Tumor Mxy', 'Healthy Mxy', 'Contrast');
      xlabel('t(s)');
      ylabel('Mxy');
      subplot(122);
      plot(evolving_t, Mzt_vxg1, 'g', evolving_t, Mzt_vxg2, 'b', evolving_t, contrast_z, 'k');
      legend('Tumor Mz', 'Healthy Mz', 'Contrast');
      xlabel('t(s)');
      ylabel('Mz');
    end
    
    function h = test_CW_dpf(obj, cw_angle, cw_phi, tr, rd_fai, t_tot, dt, rotating_angle, wf);
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      obj.init_positions();
      obj.init_magnet();
      obj.update_MeanM();
      obj.VGMean();
      Mxyt_vxg1 = zeros(1, Nt);
      Mxyt_vxg2 = zeros(1, Nt);
      Mzt_vxg1 = zeros(1, Nt);
      Mzt_vxg2 = zeros(1, Nt);
      for it = 1: Nt
        mxm_1 = obj.Sx_1; mym_1 = obj.Sy_1;
        mxm_2 = obj.Sx_2; mym_2 = obj.Sy_2;
        mxm = mxm_1 * obj.r12 + mym_2 * (1 - obj.r12);
        mym = mym_1 * obj.r12 + mym_2 * (1 - obj.r12);
        obj.RW(dt);
        obj.update_bold();
        obj.update_dpf(wf);
        obj.update_spdf();
        
        obj.CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
        obj.update_MeanM();
        obj.redeem_df();
        obj.VGMean();
        obj.flip_angle(rotating_angle);
        Mxyt_vxg1(it) = obj.Sxy_1;
        Mxyt_vxg2(it) = obj.Sxy_2;
        Mzt_vxg1(it) = obj.Sz_1;
        Mzt_vxg2(it) = obj.Sz_2;
      end
      contrast_z = abs(Mzt_vxg1 - Mzt_vxg2);
      contrast_xy = abs(Mxyt_vxg1 - Mxyt_vxg2);
      h = figure;
      title('Signal Evolving with time under CW MRI, with Bold only.');
      subplot(121);
      plot(evolving_t, Mxyt_vxg1, 'g', evolving_t, Mxyt_vxg2, 'b', evolving_t, contrast_xy, 'k');
      legend('Tumor Mxy', 'Healthy Mxy', 'Contrast');
      xlabel('t(s)');
      ylabel('Mxy');
      subplot(122);
      plot(evolving_t, Mzt_vxg1, 'g', evolving_t, Mzt_vxg2, 'b', evolving_t, contrast_z, 'k');
      legend('Tumor Mz', 'Healthy Mz', 'Contrast');
      xlabel('t(s)');
      ylabel('Mz');
    end
    
    function SE_process(obj, N_te, t_tot, dt, rotating_angle)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      obj.init_magnet();
      obj.init_positions();
      obj.update_MeanM();
      obj.SE(1);
      
      for it = 1: Nt
        obj.RW(dt);
        if it == N_te
          obj.SE(2);
        end
        obj.update_bold();
        obj.update_spdf();
        
        obj.FP_evolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
        obj.VGMean();
        obj.flip_angle(rotating_angle);
      end
    end
    
    function h = test_SEprocess(obj, N_te, t_tot, dt, rotating_angle)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      Mxyt_vxg1 = zeros(1, Nt);
      Mxyt_vxg2 = zeros(1, Nt);
      
      obj.init_magnet();
      obj.init_positions();
      obj.update_MeanM();
      obj.SE(1);
      
      for it = 1: Nt
        obj.RW(dt);
        if it == N_te
          obj.SE(2);
        end
        obj.update_bold();
        obj.update_spdf();
        
        obj.FP_evolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
        obj.VGMean();
        obj.flip_angle(rotating_angle);
        Mxyt_vxg1(it) = obj.Sxy_1;
        Mxyt_vxg2(it) = obj.Sxy_2;
      end
      contrast_xy = abs(Mxyt_vxg1 - Mxyt_vxg2);
      h = figure;
      title('Signal Evolving with time under SE MRI, with Bold only.');
      plot(evolving_t, Mxyt_vxg1, 'g', evolving_t, Mxyt_vxg2, 'b', evolving_t, contrast_xy, 'k');
      legend('Tumor Mxy', 'Healthy Mxy', 'Contrast');
      xlabel('t(s)');
    end
    
    function h = test_SE_dpf(obj, N_te, t_tot, dt, rotating_angle, wf)
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      Mxyt_vxg1 = zeros(1, Nt);
      Mxyt_vxg2 = zeros(1, Nt);
      
      obj.init_magnet();
      obj.init_positions();
      obj.update_MeanM();
      obj.SE(1);
      
      for it = 1: Nt
        obj.RW(dt);
        if it == N_te
          obj.SE(2);
        end
        obj.update_bold();
        obj.update_dpf(wf);
        obj.update_spdf();
        
        obj.FP_evolve(dt);
        obj.update_MeanM();
        obj.redeem_df();
        obj.VGMean();
        obj.flip_angle(rotating_angle);
        Mxyt_vxg1(it) = obj.Sxy_1;
        Mxyt_vxg2(it) = obj.Sxy_2;
      end
      contrast_xy = abs(Mxyt_vxg1 - Mxyt_vxg2);
      h = figure;
      title('Signal Evolving with time under SE MRI, with Bold only.');
      plot(evolving_t, Mxyt_vxg1, 'g', evolving_t, Mxyt_vxg2, 'b', evolving_t, contrast_xy, 'k');
      legend('Tumor Mxy', 'Healthy Mxy', 'Contrast');
      xlabel('t(s)');
    end
    
    function update_bold(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).update_bold(obj.Y_vxg1);
        obj.voxels_2(ivx).update_bold(obj.Y_vxg2);
      end
    end
    
    function update_dpf(obj, wf)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).update_dpf(wf);
        obj.voxels_2(ivx).update_dpf(wf);
      end
    end
    
    function update_spdf(obj)
      for ivx = 1: obj.N_voxel
        obj.voxels_1(ivx).update_spdf();
        obj.voxels_2(ivx).update_spdf();
      end
    end
    
    function redeem_df(obj)
      for ivx = 1: obj.N_voxel
        obj.voxels_1(ivx).redeem_df();
        obj.voxels_2(ivx).redeem_df();
      end
    end
    
    function update_VA(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels_1(ivx).update_VA();
        obj.voxels_2(ivx).update_VA();
      end
    end
    
    function h = plot_df(obj, plane_z)
      h = figure;
      subplot(121);
      obj.voxels_1(1).plot_df(plane_z, obj.Y_vxg1);
      subplot(122);
      obj.voxels_2(1).plot_df(plane_z, obj.Y_vxg2);
    end
    
    function flip_angle(obj, rotating_angle)
      for ivx = 1: obj.N_voxel
        if rotating_angle(1)
          obj.voxels_1(ivx).flip_angle();
        end
        if rotating_angle(2)
          obj.voxels_2(ivx).flip_angle();
        end
      end
    end
  end%methods end
  
end%class end
    
    