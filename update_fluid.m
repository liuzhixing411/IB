function Fluid=update_fluid(u_pred,v_pred,alpha_new,f_x,f_y,Fluid,grid)
% 更新流体速度和压力
% 步骤：
% 1. u_new = u_pred + f*dt (施加相互作用力)
% 2. 求解压力泊松方程
% 3. 速度修正（在泊松方程中已经完成）

dt=grid.dt;
ghostnum=grid.ghostnum;
h=grid.h;
u_old=Fluid.u;
v_old=Fluid.v;
rhol=Fluid.rhol;
rhog=Fluid.rhog;
Ny=grid.Ny;
Nx=grid.Nx;
endy=grid.endy;
endx=grid.endx;


%% 步骤1：应用相互作用力到预测速度
u_new=u_pred+f_x*dt;
v_new=v_pred+f_y*dt;

%% 边界条件（无滑移边界）
% left ghosts for u
for j = 1:(Ny+1)
    for i = 1:ghostnum
        u_new(j,i) = - u_new(j, 2*ghostnum + 2 - i);
    end
end

% right ghosts for u
for j = 1:(Ny+1)
    for i = endx+1:(Nx+1)
        u_new(j,i) = - u_new(j, 2*endx + 2 - i);
    end
end

% top ghosts for u
for j = 1:ghostnum
    for i = 1:(Nx+1)
        u_new(j,i) = - u_new(2*ghostnum + 1 - j, i);
    end
end

% bottom ghosts for u
for j = endy+1:(Ny+1)
    for i = 1:(Nx+1)
        u_new(j,i) = - u_new(2*endy + 1 - j, i);
    end
end

% 壁面u速度设为0
u_new(:, ghostnum+1) = 0;
u_new(:, endx+1) = 0;

% left ghosts for v
for j = 1:(Ny+1)
    for i = 1:ghostnum
        v_new(j,i) = - v_new(j, 2*ghostnum + 1 - i);
    end
end

% right ghosts for v
for j=1:Ny+1
    for i=endx+1:Nx+1
        v_new(j,i) = - v_new(j,2*endx+1-i);
    end
end

% top ghosts for v
for j = 1:ghostnum
    for i = 1:(Nx+1)
        v_new(j,i) = - v_new(2*ghostnum + 2 - j, i);
    end
end

% bottom ghosts for v
for j = endy+1:(Ny+1)
    for i = 1:(Nx+1)
        v_new(j,i) = - v_new(2*endy + 2 - j, i);
    end
end

% 壁面v速度设为0
v_new(ghostnum+1, :) = 0;
v_new(endy+1, :) = 0;

%% 步骤2：求解压力泊松方程
% 计算混合密度
rho_cc = rhol .* alpha_new + rhog .* (1 - alpha_new);

% 调用泊松求解器
% 注意：Poisson函数内部会进行速度修正
% 但由于我们已经在prediction中包含了n时刻的压力梯度
% 这里求解的是n+1时刻的压力修正
p_new=Poisson(rho_cc,u_new,v_new,u_old,v_old,h,h,dt,ghostnum,f_x,f_y);

%% 步骤3：压力修正（如果需要）
% 在当前实现中，压力修正已经在下一个时间步的prediction中完成
% 这里直接使用u_new和v_new

%% 赋值到Fluid结构体
Fluid.u=u_new;
Fluid.v=v_new;
Fluid.alpha=alpha_new;
Fluid.p=p_new;

end