function Particle=update_particle(f_x,f_y,Particle,grid)
% 更新粒子的平移和旋转运动
% 根据论文 Eq.(28)-(33)

g=9.8;
dt=grid.dt;
start=grid.start;
endy=grid.endy;
endx=grid.endx;
h=grid.h;
ghostnum = grid.ghostnum;

% 粒子参数
u=Particle.u;
v=Particle.v;
omega=Particle.omega;
m=Particle.m;
x_c=Particle.x_c;
y_c=Particle.y_c;

% 计算转动惯量 (2D圆盘)
I = Particle.I;

% 计算流体作用在粒子上的力和力矩
% 根据论文 Eq.(28): F = -∫f dx (注意负号!)
% 但是f已经是流体对粒子的作用力，所以粒子受到的力应该是-F
% 实际上，在IBM中，f是施加在流体上使其满足边界条件的力
% 根据牛顿第三定律，粒子受到的力是-f

% 重新计算：
Fx = 0;
Fy = 0;
Torque = 0;

for j=start:endy+1
    for i=start:endx+1
        % u-face: 贡献f_x和相应的力矩
        if i >= start && i <= endx+1
            
            y_u = (j - ghostnum -0.5) * h;
            Fx = Fx - f_x(j,i) * h^2;
            % 力矩：r × F，2D中T_z = r_x * F_y - r_y * F_x
            % 但f_x在u-face上，我们需要v的分量来计算力矩
            % 简化：只使用y方向的力矩贡献
            Torque = Torque + (y_u - y_c) * f_x(j,i) * h^2;
        end
        
        % v-face: 贡献f_y和相应的力矩  
        if j >= start && j <= endy+1
            x_v = (i - ghostnum -0.5) * h;
            
            Fy = Fy - f_y(j,i) * h^2;
            % 力矩：r × F
            Torque = Torque - (x_v - x_c) * f_y(j,i) * h^2;
        end
    end
end

% 更新速度和角速度 - 根据论文 Eq.(30) 和 Eq.(32)
% dU/dt = F/M + g
% dω/dt = T/I

u_new = u + (Fx/m) * dt;
v_new = v + (Fy/m - g) * dt;  % 重力向下
omega_new = omega + (Torque/I) * dt;


% 更新位置和角度 - 根据论文 Eq.(31) 和 Eq.(33)
% dX/dt = U (使用中点法或显式Euler)
x_c_new = x_c + 0.5*(u + u_new) * dt;
y_c_new = y_c + 0.5*(v + v_new) * dt;

% 更新粒子结构体
Particle.u = u_new;
Particle.v = v_new;
Particle.omega = omega_new;
Particle.x_c = x_c_new;
Particle.y_c = y_c_new;

end