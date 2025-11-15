function Fluid=update_fluid(u_pred,v_pred,alpha_new,f_x,f_y,Fluid,grid)

dt=grid.dt;

%% 步骤1：预测速度

u_star = u_pred + f_x*dt;
v_star = v_pred + f_y*dt;

%% 步骤2：求解压力泊松方程并修正速度

Fluid.alpha = alpha_new; 
[u_new, v_new, p_new] = Poisson(Fluid, grid, u_star, v_star);

%% 步骤3：赋值到Fluid结构体
Fluid.u = u_new;
Fluid.v = v_new;
Fluid.alpha = alpha_new;
Fluid.p = p_new; 
end
