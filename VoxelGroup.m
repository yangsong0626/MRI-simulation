%Voxel group is a set of voxels with different properties

classdef VoxelGroup<handle
  properties
    %Voxels
    voxels
    N_voxel
    
    %Voxel posistion
    voxel_x
    voxel_y
    voxel_z
    
    voxel_angle
    
    %Magnetization signal
    Sx
    Sy
    Sz
    Sxy
  end%properties end
  
  methods
    function obj = VoxelGroup(N_voxel, x, num_spin, v, T1, T2, diameter, rsize, Diffusion, theta)
      obj.voxels = VXL.empty(0, N_voxel);
      for ivx = 1: N_voxel
          obj.voxels(ivx) = VXL(num_spin, v, T1, T2, diameter, rsize, Diffusion, theta(ivx));
      end
      obj.N_voxel = N_voxel;
      
      obj.voxel_x = x; 
    end
    
    function VGMean(obj)
      obj.Sx = 0;
      obj.Sy = 0;
      obj.Sz = 0;
      obj.Sxy = 0;
      for ivx = 1: obj.N_voxel
        obj.Sx = obj.Sx + obj.voxels(ivx).Mean_Mx / obj.N_voxel;
        obj.Sy = obj.Sy + obj.voxels(ivx).Mean_My / obj.N_voxel;
        obj.Sxy = obj.Sxy + obj.voxels(ivx).Mxy / obj.N_voxel;
        obj.Sz = obj.Sz + obj.voxels(ivx).Mean_Mz / obj.N_voxel;
      end
    end
    
    function init_magnet(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).init_magnet();
      end
    end
    
    function init_positions(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).init_positions();
      end
    end
    
    function RW(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).RW();
      end
    end
    
    function update_magnet(obj, dM)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).update_magnet(dM);
      end
    end
    
    function update_MeanM(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).update_MeanM();
      end
    end
    
    function update_df(obj, df)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).update_df(df);
      end
    end
    
    function SE(obj, pls_ord)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).SE(pls_ord);
      end
    end
    
    function PNNP(obj, pls_ord, pls_angle)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).PNNP(pls_ord, pls_angle);
      end
    end
    
    function IR(obj, pls_ord)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).IR(pls_ord);
      end
    end
    
    function SS(obj, pls_angle)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).SS(pls_angle);
      end
    end
    
    function RD_evolve(obj, mxm, mym, tr, fai, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).RD_evolve(mxm, mym, tr, fai, dt);
      end
    end
    
    function FP_evolve(obj, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).FP_evolve(dt);
      end
    end
    
    function CW_evolve(obj, cw_angle, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).CW_evolve(cw_angle, dt);
      end
    end
    
    function CWRD_evolve(obj, cw_angle, mxm, mym, tr, rd_fai, dt)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).CWRD_evolve(cw_angle, mxm, mym, tr, rd_fai, dt);
      end
    end
    
    function CWRD_RMevolve(obj, cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt)
      for ivx = 1: obj.N_voxel
        obj.voxels(ivx).CWRD_RMevolve(cw_angle, cw_phi, mxm, mym, tr, rd_fai, dt);
      end
    end
    
    function CWRD_RMprocess(obj, cw_angle, cw_phi, mxm, mym, tr, rd_fai, t_tot, dt, Y_vx, rotating_angle);
      for ivx = 1: obj.N_voxel
        obj.voxels(ivx).CWRD_RMprocess(cw_angle, cw_phi, mxm, mym, tr, rd_fai, t_tot, dt, Y_vx, rotating_angle);
      end
    end
    
    function SE_process(obj, N_te, Y_vx, t_tot, dt, rotating_angle)
      for ivx = 1: obj.N_voxel
        obj.voxels(ivx).SE_process(N_te, Y_vx, t_tot, dt, rotating_angle);
      end
    end
    
    function update_bold(obj, Y)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).update_bold(Y);
      end
    end
    
    function update_dpf(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).update_dpf();
      end
    end
    
    function update_spdf(obj)
      for ivx = 1: obj.N_voxel
        obj.voxels(ivx).update_spdf();
      end
    end
    
    function redeem_df(obj)
      for ivx = 1: obj.N_voxel
        obj.voxels(ivx).redeem_df();
      end
    end
    
    function update_VA(obj)
      for ivx = 1 : obj.N_voxel
        obj.voxels(ivx).update_VA();
      end
    end
    
    function flip_angle(obj)
      for ivx = 1: obj.N_voxel
        obj.voxels(ivx).flip_angle();
      end
    end
  end%methods end
  
end%class end
    
    