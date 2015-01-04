%class to simulate double-spin system
classdef DoubleSpin<handle
  properties
    Spin1
    Spin2
    tr
    RD_phase
    w1
    CW_phase
    
    %relaxation constant
    T1 = Inf;
    T2 = Inf;
    
    %ratio of contribution to feedback
    r12;
  end
  
  methods
    function obj = DoubleSpin(df1, df2, tr, RD_phase, w1, CW_phase, r12)
      obj.Spin1.M = [0; 0; 1];
      obj.Spin2.M = [0; 0; 1];
      obj.Spin1.df = df1;
      obj.Spin2.df = df2;
      
      obj.tr = tr;
      obj.RD_phase = RD_phase;
      
      obj.w1 = w1;
      obj.CW_phase = CW_phase;
      
      obj.r12 = r12;
    end
    function init_M(obj)
      obj.Spin1.M = [0; 0; 1];
      obj.Spin2.M = [0; 0; 1];
    end
    
    function CWRD_RMevolve(obj, dt)
      mxm = obj.Spin1.M(1) * obj.r12 + obj.Spin2.M(1) * (1 - obj.r12);
      mym = obj.Spin1.M(2) * obj.r12 + obj.Spin2.M(2) * (1 - obj.r12);
      
      rotmat = rdrot(mxm, mym, obj.RD_phase, obj.w1, obj.CW_phase, obj.Spin1.df, obj.tr, dt);
      obj.Spin1.M = rotmat * obj.Spin1.M;
      
      rotmat = rdrot(mxm, mym, obj.RD_phase, obj.w1, obj.CW_phase, obj.Spin2.df, obj.tr, dt);
      obj.Spin2.M = rotmat * obj.Spin2.M;
    end
    
    function CWRD_process(obj, t_tot, dt)
      obj.init_M();
      evolving_t = 0: dt: t_tot;
      Nt = length(evolving_t);
      
      Mt_1 = zeros(3, Nt + 1);
      Mt_2 = zeros(3, Nt + 1);
      
      Mt_1(:, 1) = obj.Spin1.M;
      Mt_2(:, 2) = obj.Spin2.M;
      
      for it = 1: Nt
        obj.CWRD_RMevolve(dt);
        Mt_1(:, it + 1) = obj.Spin1.M;
        Mt_2(:, it + 1) = obj.Spin2.M;
      end
      
      [x y z] = sphere(128);
      h1 = surf(x, y, z); 
      set(h1, 'FaceAlpha', 0, 'MeshStyle', 'column', 'EdgeAlpha', 0.1)
      hold on;
      plot3(Mt_1(1, :), Mt_1(2, :), Mt_1(3, :), 'b', Mt_2(1, :), Mt_2(2, :), Mt_2(3, :), 'g', 'LineWidth', 1.5);
      xlabel('Mx');
      ylabel('My');
      zlabel('Mz');
      title('Magnetization evolve');
      hold off;
    end
    
    function F = CWphase(obj, t_tot, dt, Nframe)
      %reset CW_phase to 0
      obj.CW_phase = 0;
      %save gif into this file
      filename = 'influence of CW phase.gif';
      %prepare the frames
      F(Nframe) = struct('cdata', [], 'colormap', []);
      dCWphase = 2 * pi / Nframe;
      h = figure;
      for iframe = 1: Nframe
        refreshdata(h);
        obj.CWRD_process(t_tot, dt);
        drawnow
        pause(0.02);
        F(iframe) = getframe(gcf);
        im = frame2im(F(iframe));
        [A, map] = rgb2ind(im, 256);
        if iframe == 1
          imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1);
        else
          imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
        end
        obj.CW_phase = obj.CW_phase + dCWphase;
      end
      movie(F);
    end
      
    function CWstrength(obj, t_tot, dt, Nframe)
      obj.w1 = 0;
        %prepare the frames
      F(1: Nframe) = struct('cdata', [], 'colormap', []);
      dw1 = 100 / Nframe;
      for iframe = 1: Nframe
        obj.CWRD_process(t_tot, dt);
        F(iframe) = getframe(gcf);
        obj.w1 = obj.w1 + dw1;
      end
      movie(F);
      
      filename = 'influence of CW strength.gif';
      for iframe = 1: Nframe
        im = frame2im(F(iframe));
        [A, map] = rgb2ind(im, 256);
        if iframe == 1
          imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1);
        else
          imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
        end
      end
      
    end
  end%method end
end%class end

      