function [u_new, v_new, p_new] = Poisson(Fluid, grid, u_star, v_star)

% 1. 完整实现变密度泊松方程 
%    LHS: ∇ ⋅ ( (1/ρ) ∇p )
%    RHS: (1/Δt) ∇ ⋅ u*
% 2. 假设单位厚度 (h_depth = 1.0)
% 3. 速度修正 u_new = u* - (Δt/ρ) ∇p

% 0. 获取参数
alpha = Fluid.alpha;
rhol = Fluid.rhol; 
rhog = Fluid.rhog; 
p_old = Fluid.p;   

dt = grid.dt;
h = grid.h;
ghostnum = grid.ghostnum;
Ny = grid.Ny;
Nx = grid.Nx;
start = grid.start;
endy = grid.endy;
endx = grid.endx;

eps_div = 1e-12;
h2 = h * h;

% --- 1. 计算密度 (假设单位厚度) ---
h_depth = 1.0; % [m]


rho_cc_3D = rhol .* alpha + rhog .* (1 - alpha);

rho_cc = rho_cc_3D * h_depth; 


rhohalfx = zeros(Ny+1, Nx+1);
for j=start:endy+1
    for i=start:endx+1
        rhohalfx(j,i) = (rho_cc(j,i) + rho_cc(j,i-1))/2;
    end
end

rhohalfy = zeros(Ny+1, Nx+1);
for j=start:endy+1
    for i=start:endx+1
        rhohalfy(j,i) = (rho_cc(j,i) + rho_cc(j-1,i))/2;
    end
end

% 2. 计算 RHS (Sval) 
%  RHS = (1/Δt) ∇ ⋅ u*
Sval = zeros(Ny, Nx);
for j=start:endy
    for i=start:endx
        % div(u*) 在 cell center (j,i)
        div_u_star = (u_star(j, i+1) - u_star(j, i))/h + ...
                     (v_star(j+1, i) - v_star(j, i))/h;
        
        Sval(j,i) = (1.0 / dt) * div_u_star;
    end
end

% --- 3. 迭代求解泊松方程 ∇ ⋅ ( (1/ρ) ∇p ) = Sval ---
p_new = p_old; % 使用上一步压力作为初始猜测
p_iter = p_old;
max_iter = 5000; 
tol = 1e-5;

for iter = 1:max_iter
    p_prev_iter = p_iter;

    % -------------------------
    % 1) 在使用 p_iter 进行更新前，先为 p_iter 设置 ghost 值（基于 u_star, v_star）
    % Left boundary (ghost i = start-1)
    for j = start:endy
        rho_face = rhohalfx(j, start);        
        uface = u_star(j, start);             
        p_iter(j, start-1) = p_iter(j, start) - h * (rho_face / dt) * (uface - 0); % u_b = 0
    end
    % Right boundary (ghost i = endx+1)
    for j = start:endy
        rho_face = rhohalfx(j, endx+1);
        uface = u_star(j, endx+1);
        p_iter(j, endx+1) = p_iter(j, endx) - h * (rho_face / dt) * (uface - 0);
    end
    % Bottom boundary (ghost j = start-1)
    for i = start:endx
        rho_face = rhohalfy(start, i);
        vface = v_star(start, i);
        p_iter(start-1, i) = p_iter(start, i) - h * (rho_face / dt) * (vface - 0);
    end
    % Top boundary (ghost j = endy+1)
    for i = start:endx
        rho_face = rhohalfy(endy+1, i);
        vface = v_star(endy+1, i);
        p_iter(endy+1, i) = p_iter(endy, i) - h * (rho_face / dt) * (vface - 0);
    end
    % -------------------------
    
    % -------------------------
    % 2) 用 p_iter（现在已含与 BC 一致的 ghost）计算 p_new（Jacobi 更新）
    for j = start:endy
        for i = start:endx
            Ae = 1.0 / (rhohalfx(j, i+1) + eps_div);
            Aw = 1.0 / (rhohalfx(j, i)   + eps_div);
            An = 1.0 / (rhohalfy(j+1, i) + eps_div);
            As = 1.0 / (rhohalfy(j, i)   + eps_div);
            Ap = Ae + Aw + An + As;
            p_new(j,i) = ( (Ae * p_iter(j,i+1) + Aw * p_iter(j,i-1)) + ...
                           (An * p_iter(j+1,i) + As * p_iter(j-1,i)) - ...
                           h2 * Sval(j,i) ) / Ap;
        end
    end
    % -------------------------

    % （可选）去掉压力平均值以消除纯 Neumann 导致的常数模
    interior = p_new(start:endy, start:endx);
    mean_p = mean(interior(:));
    p_new(start:endy, start:endx) = interior - mean_p;

    % 更新 p_iter，准备下一次迭代
    p_iter = p_new;

    % 计算收敛残差
    res_norm = norm(p_iter(start:endy, start:endx) - p_prev_iter(start:endy, start:endx), 'fro');
    p_norm = norm(p_iter(start:endy, start:endx), 'fro') + eps_div;
    res = res_norm / p_norm;
    if res < tol
        break;
    end
end
% 最终 p_new = p_iter;
p_new = p_iter;

% if iter == max_iter
%     fprintf('Warning: Poisson did not converge. Res = %e\n', res);
% end

%  4. 速度修正
% u_new = u_star - (Δt/ρ_face) * ∇p
u_new = zeros(Ny+1, Nx+1);
v_new = zeros(Ny+1, Nx+1);

% u-velocity
for j=start:endy+1 
    for i=start+1:endx+1
        grad_p_x = (p_new(j,i) - p_new(j,i-1)) / h;
        % 使用 u-face 上的 2D 密度 rhohalfx
        u_new(j,i) = u_star(j,i) - (dt / (rhohalfx(j,i) + eps_div)) * grad_p_x;
    end
end

% v-velocity
for j=start+1:endy+1 
    for i=start:endx+1
        grad_p_y = (p_new(j,i) - p_new(j-1,i)) / h;
        % 使用 v-face 上的 2D 密度 rhohalfy
        v_new(j,i) = v_star(j,i) - (dt / (rhohalfy(j,i) + eps_div)) * grad_p_y;
    end
end

% 5. 应用速度边界条件 (No-Slip Walls)

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

end
