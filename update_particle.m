function Particle=update_particle(f_x,f_y,Particle,Fluid,grid)
% 更新粒子的平移和旋转运动

% m*dU/dt = ∮(τ ⋅ n_p) dA + (ρ_p - ρ_f) * V_p * g

%  0. 获取参数

rhol = Fluid.rhol; 
rhog = Fluid.rhog; 
alpha = Fluid.alpha;
g_val = 9.8;
dt = grid.dt;
h = grid.h;
ghostnum = grid.ghostnum;

% 粒子参数
rhos= Particle.rho;

Vp_area = Particle.V;

m = Particle.m;

I = Particle.I;

x_c = Particle.x_c;
y_c = Particle.y_c;

%相互作用力

rho_cc = rhol .* alpha + rhog .* (1 - alpha);
% 网格单元的体积 dV
dV = h * h ; 

Fx = 0;
Fy = 0;
Torque = 0;

for j=grid.start:grid.endy+1
    for i=grid.start:grid.endx+1
        
        % u-face (j, i)
        if i >= grid.start && i <= grid.endx+1
            
            rho_face_u_3D = (rho_cc(j,i) + rho_cc(j,i-1))/2;
            y_u = (j - ghostnum - 0.5) * h;            
            Fx_contrib = - (rho_face_u_3D * f_x(j,i) * dV);
            Fx = Fx + Fx_contrib;
            Torque = Torque - (y_u - y_c) * Fx_contrib; 
        end
        
        % v-face (j, i)
        if j >= grid.start && j <= grid.endy+1
            
            rho_face_v_3D = (rho_cc(j,i) + rho_cc(j-1,i))/2;
            x_v = (i - ghostnum - 0.5) * h;
            Fy_contrib = - (rho_face_v_3D * f_y(j,i) * dV);
            Fy = Fy + Fy_contrib;
            Torque = Torque + (x_v - x_c) * Fy_contrib;
        end
    end
end

% 3. 计算净浮力
% (ρ_p - ρ_f) * V_p * g 

% 粒子体积 V_p
Vp_volume = Vp_area; 
rho_f_local_3D = Fluid.rhol; 
F_buoyancy_y = (rhos - rho_f_local_3D) * Vp_volume * (-g_val);
F_buoyancy_x = 0; % 假设重力只在y方向

% 4. 更新运动
% m*dU/dt = F_hydro + F_buoyancy

% 粒子质量 m
m_particle = m ; % [kg]
% 粒子转动惯量 I
I_particle = I ; % [kg*m^2]
% (检查: I_init = 0.5*m*r^2 = 0.5*(m_particle/h_depth)*r^2 = I_particle/(h_depth*2) * r^2 ... 


% 总力
F_total_x = Fx + F_buoyancy_x;
F_total_y = Fy + F_buoyancy_y;

% 更新速度 (一步欧拉法)
u_new = Particle.u + (F_total_x / m_particle) * dt;
v_new = Particle.v + (F_total_y / m_particle) * dt;
omega_new = Particle.omega + (Torque / I_particle) * dt;

% 更新位置
x_c_new = x_c + u_new * dt;
y_c_new = y_c + v_new * dt;

%  5. 赋值 
Particle.u = u_new;
Particle.v = v_new;
Particle.omega = omega_new;
Particle.x_c = x_c_new;
Particle.y_c = y_c_new;
end
