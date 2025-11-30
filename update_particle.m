function Particle = update_particle(f_x, f_y, Particle, Fluid, grid, Wall)
% UPDATE_PARTICLE
% 更新粒子运动，包含：
% 1. 流体动力 (IBM 宏观力)
% 2. 润滑力修正 (Lubrication Correction, Paper Sec. III)
% 3. 软球碰撞模型 (Soft-sphere Collision, Paper Sec. II)
% 4. 子步时间积分 (Sub-stepping, Paper Sec. IV)
%
% 参考文献: Costa et al., "Collision model for fully-resolved simulations...", 2015.

%% 1. 物理参数定义 (根据论文 Table II 和 Table IV 设置)
% -----------------------------------------------------------
N_col   = 8;       % 碰撞时间倍数 T_n = N_col * dt (Paper Eq. 18)
e_n_d   = 0.97;    % 干法向恢复系数 (restitution coefficient)
e_t_d   = 0.39;    % 干切向恢复系数 (tangential restitution)
mu_c    = 0.15;    % 摩擦系数 (Coulomb friction)
K_gyr   = 1/2;     

% 润滑力参数 (Paper Sec. III)
eps_dx  = 0.05;    % IBM 失效的间隙阈值 (epsilon_Delta_x)
eps_sigma = 0.001; % 粗糙度截断阈值 (epsilon_sigma)
mu_fluid = 0.001;  % 流体动力粘度 (Pa.s), 需根据实际流体设置，此处假设水

% 子步积分参数 (Paper Sec. IV)
N_sub   = 10;      % 将一个流体时间步 dt 分为 N_sub 个子步进行粒子积分
dt_sub  = grid.dt / N_sub;

%% 2. 准备工作
Np = length(Particle.r);
if ~isfield(Particle, 'delta_t') || isempty(Particle.delta_t)
    Particle.delta_t = cell(Np, Np);      % 粒子间切向位移记忆
    Particle.n_prev  = cell(Np, Np);      % 上一时刻法向
    Particle.delta_t_wall = cell(Np, 4);  % 墙面切向位移记忆
end

% 预计算流体网格属性
dV = grid.h * grid.h;
rho_mix = Fluid.rhol .* Fluid.alpha + Fluid.rhog .* (1 - Fluid.alpha);
g_vec = [0; -9.81]; % 重力加速度

%% 3. 计算流体宏观力 (F_hydro)
% 注意：F_hydro 在整个 dt 时间步内假设为常数 (Time splitting)
F_hydro_x = zeros(Np, 1);
F_hydro_y = zeros(Np, 1);
T_hydro   = zeros(Np, 1);

for p = 1:Np
    x_c = Particle.x_c(p);
    y_c = Particle.y_c(p);

    % 简单的离散力求和 (根据传入的 f_x, f_y)
    % 实际应用中通常使用插值核函数 (Delta function) 将力投影回粒子
    % 此处保留原代码的简化逻辑，遍历包围盒
    start_i = max(1, floor((x_c - Particle.r(p))/grid.h));
    end_i   = min(grid.endx+1, ceil((x_c + Particle.r(p))/grid.h));
    start_j = max(1, floor((y_c - Particle.r(p))/grid.h));
    end_j   = min(grid.endy+1, ceil((y_c + Particle.r(p))/grid.h));

    fx_sum = 0; fy_sum = 0; t_sum = 0;

    for j = start_j:end_j
        for i = start_i:end_i
            % 简单的最近邻搜索或核函数积分范围
            % 此处仅作示意，沿用原逻辑框架
            if i > 1 && i <= size(f_x,2) && j <= size(f_x,1)
                rho_u = (rho_mix(j,i) + rho_mix(j,i-1))/2;
                force_val = -(rho_u * f_x(j,i) * dV);
                fx_sum = fx_sum + force_val;
                % 力矩臂 y
                dist_y = (j - grid.ghostnum - 0.5)*grid.h - y_c;
                t_sum = t_sum - dist_y * force_val;
            end

            if j > 1 && j <= size(f_y,1) && i <= size(f_y,2)
                rho_v = (rho_mix(j,i) + rho_mix(j-1,i))/2;
                force_val = -(rho_v * f_y(j,i) * dV);
                fy_sum = fy_sum + force_val;
                % 力矩臂 x
                dist_x = (i - grid.ghostnum - 0.5)*grid.h - x_c;
                t_sum = t_sum + dist_x * force_val;
            end
        end
    end
    F_hydro_x(p) = fx_sum;
    F_hydro_y(p) = fy_sum;
    T_hydro(p)   = t_sum;
end

%% 4. 子步积分循环 (Sub-stepping Loop)
% 论文指出，碰撞力和润滑力需要在更小的时间尺度上积分

for step = 1:N_sub
    for p = 1:Np
        % 提取当前子步状态
        m = Particle.m(p);
        I = Particle.I(p);
        r_p = Particle.r(p);
        pos_p = [Particle.x_c(p); Particle.y_c(p)];
        vel_p = [Particle.u(p); Particle.v(p)];
        omega_p = Particle.omega(p);

        % 初始化碰撞力与润滑力
        F_cont_x = 0; F_cont_y = 0;
        T_cont = 0;

        % --- A. 粒子-粒子 交互 (Collision & Lubrication) ---
        for q = 1:Np
            if p == q, continue; end

            pos_q = [Particle.x_c(q); Particle.y_c(q)];
            vel_q = [Particle.u(q); Particle.v(q)];
            omega_q = Particle.omega(q);
            r_q = Particle.r(q);
            m_q = Particle.m(q);

            rel_pos = pos_q - pos_p;
            dist = norm(rel_pos);
            if dist == 0, continue; end
            n_vec = rel_pos / dist; % 法向量 p -> q

            gap = dist - (r_p + r_q);
            eps = gap / r_p; % 归一化间隙 (假设等径或以 p 为主)

            % 接触点相对速度 (Eq. 11)
            % u_cp = u + omega x R
            v_surf_p = vel_p + cross_product_2d(omega_p, r_p * n_vec);
            v_surf_q = vel_q + cross_product_2d(omega_q, -r_q * n_vec);
            u_ij = v_surf_p - v_surf_q;

            u_n_val = dot(u_ij, n_vec);
            u_n_vec = u_n_val * n_vec;
            u_t_vec = u_ij - u_n_vec;

            % 1. 润滑力计算 (Lubrication Force, Section III)
            % 仅在靠近且未接触，且相互接近时计算
            if gap > 0 && gap < eps_dx * grid.h && u_n_val > 0
                % 使用论文 Eq. 32 (粒子-粒子) 或 Eq. 33 (粒子-墙)
                % 此处简化为两球模型
                lambda_eps = lubrication_factor_pp(eps, eps_sigma, eps_dx * grid.h / r_p);
                % F_lub = -6 * pi * mu * R * u_n * (correction)
                F_lub_val = -6 * pi * mu_fluid * r_p * u_n_val * lambda_eps;
                F_lub = F_lub_val * n_vec;

                F_cont_x = F_cont_x + F_lub(1);
                F_cont_y = F_cont_y + F_lub(2);
                % 润滑力沿法向，无力矩
            end

            % 2. 接触力计算 (Collision Model, Section II)
            overlap = -gap;
            if overlap > 0
                % 等效质量 (Eq. 19)
                m_eff = (m * m_q) / (m + m_q);

                % 计算刚度和阻尼系数 (Eq. 18, 25)
                % 碰撞时间 T_n = N_col * dt (注意是大时间步 dt)
                [kn, eta_n] = calc_coeffs(m_eff, e_n_d, N_col * grid.dt);

                % 法向力 (Eq. 15)
                F_n = -kn * overlap * n_vec - eta_n * u_n_vec;
                % 确保只有挤压(排斥)力，没有"拉"力（虽然 overlap>0 隐含挤压，但阻尼项可能导致反向）
                % 通常软球模型直接应用公式即可。

                % 切向力处理 (Eq. 22-28)
                m_eff_t = m_eff / (1 + 1/K_gyr);
                [kt, eta_t] = calc_coeffs(m_eff_t, e_t_d, N_col * grid.dt);

                % 读取/初始化切向历史
                % 正确获取有序标量索引（保证 idx1 < idx2 且为标量）
                idxs = sort([p, q]);
                idx1 = idxs(1);
                idx2 = idxs(2);

                if isempty(Particle.delta_t{idx1, idx2})
                    Particle.delta_t{idx1, idx2} = [0;0];
                    Particle.n_prev{idx1, idx2} = n_vec;
                end

                % 旋转上一时刻的切向位移 (Eq. 27)
                n_prev = Particle.n_prev{idx1, idx2};
                delta_t_old = Particle.delta_t{idx1, idx2};
                % 2D 旋转: 投影到旧切向 -> 映射到新切向
                t_vec_curr = [-n_vec(2); n_vec(1)];
                t_vec_prev = [-n_prev(2); n_prev(1)];
                delta_scalar = dot(delta_t_old, t_vec_prev);
                delta_t_rot = delta_scalar * t_vec_curr;

                % 积分切向位移
                delta_t_star = delta_t_rot + u_t_vec * dt_sub;

                % 试探力
                F_t_trial = -kt * delta_t_star - eta_t * u_t_vec;
                Fn_mag = norm(F_n);
                Ft_mag = norm(F_t_trial);

                % Coulomb 饱和判断 (Eq. 28)
                if Ft_mag <= mu_c * Fn_mag
                    F_t = F_t_trial;
                    Particle.delta_t{idx1, idx2} = delta_t_star; % Stick
                else
                    % Sliding
                    t_dir = F_t_trial / (Ft_mag + 1e-12);
                    F_t = mu_c * Fn_mag * t_dir; % 方向与试探力相同(即阻碍运动)
                    % 修正 delta_t 以保持一致性
                    Particle.delta_t{idx1, idx2} = (-F_t - eta_t * u_t_vec) / kt;
                end
                Particle.n_prev{idx1, idx2} = n_vec; % 更新法向

                % 累加力
                F_total_ij = F_n + F_t;
                F_cont_x = F_cont_x + F_total_ij(1);
                F_cont_y = F_cont_y + F_total_ij(2);

                % 累加力矩 (r_p x F_t)
                % 力作用点在表面: arm = r_p * n_vec
                T_cont = T_cont + (r_p * n_vec(1) * F_total_ij(2) - r_p * n_vec(2) * F_total_ij(1));
            else
                % 正确获取有序标量索引
                idxs = sort([p, q]);
                idx1 = idxs(1);
                idx2 = idxs(2);

                % 安全清除：先判断是否存在再清除（避免对空索引/非标量再次报错）
                if ~isempty(Particle.delta_t) && numel(Particle.delta_t) >= (idx1-1)*Np + idx2
                    Particle.delta_t{idx1, idx2} = [];
                end
                if ~isempty(Particle.n_prev) && numel(Particle.n_prev) >= (idx1-1)*Np + idx2
                    Particle.n_prev{idx1, idx2} = [];
                end


            end
        end

        % --- B. 粒子-墙 交互 ---
        % Walls: left, right, down, up
        walls = [Wall.left, Wall.right, Wall.down, Wall.up];
        normals = {[1;0], [-1;0], [0;1], [0;-1]}; % 指向流体内部

        for w = 1:4
            n_w = normals{w};
            % 计算到墙距离
            if w<=2, dist_w = abs(pos_p(1) - walls(w));
            else,    dist_w = abs(pos_p(2) - walls(w)); end

            gap = dist_w - r_p;
            eps = gap / r_p;

            % 墙面相对速度 (假设墙静止)
            v_surf_p = vel_p + cross_product_2d(omega_p, r_p * (-n_w)); % 接触点法向为 -n_w
            u_pw = v_surf_p; % - 0
            u_n_val = dot(u_pw, -n_w); % 逼近速度 (>0 为靠近)

            % 1. 墙面润滑力 (Paper Eq. 33)
            if gap > 0 && gap < eps_dx * grid.h && u_n_val > 0
                lambda_eps = lubrication_factor_pw(eps, eps_sigma, eps_dx * grid.h / r_p);
                F_lub_val = -6 * pi * mu_fluid * r_p * u_n_val * lambda_eps;
                % 润滑力方向沿法向 (-n_w 为推开粒子的方向)
                % 注意: 这里 n_w 是指向域内的。排斥力应沿 n_w 方向。
                F_lub = F_lub_val * n_w;
                F_cont_x = F_cont_x + F_lub(1);
                F_cont_y = F_cont_y + F_lub(2);
            end

            % 2. 墙面接触力
            overlap = -gap;
            if overlap > 0
                m_eff = m; % 墙质量无穷大
                [kn, eta_n] = calc_coeffs(m_eff, e_n_d, N_col * grid.dt);

                u_n_vec = dot(u_pw, n_w) * n_w;
                u_t_vec = u_pw - u_n_vec;

                F_n = kn * overlap * n_w - eta_n * u_n_vec; % 弹力指向域内(n_w)

                % 切向
                m_eff_t = m_eff / (1 + 1/K_gyr);
                [kt, eta_t] = calc_coeffs(m_eff_t, e_t_d, N_col * grid.dt);

                if isempty(Particle.delta_t_wall{p, w})
                    Particle.delta_t_wall{p, w} = [0;0];
                end
                delta_t_star = Particle.delta_t_wall{p, w} + u_t_vec * dt_sub;

                F_t_trial = -kt * delta_t_star - eta_t * u_t_vec;
                if norm(F_t_trial) <= mu_c * norm(F_n)
                    F_t = F_t_trial;
                    Particle.delta_t_wall{p, w} = delta_t_star;
                else
                    t_dir = F_t_trial / (norm(F_t_trial)+1e-12);
                    F_t = mu_c * norm(F_n) * t_dir;
                    Particle.delta_t_wall{p, w} = (-F_t - eta_t * u_t_vec) / kt;
                end

                F_total = F_n + F_t;
                F_cont_x = F_cont_x + F_total(1);
                F_cont_y = F_cont_y + F_total(2);

                % 力矩 (接触点矢量 r_c = -r_p * n_w)
                vec_arm = -r_p * n_w;
                T_cont = T_cont + (vec_arm(1)*F_total(2) - vec_arm(2)*F_total(1));
            else
                Particle.delta_t_wall{p, w} = [];
            end
        end

        % --- C. 子步更新 (积分) ---
        % 合力 = F_hydro + F_gravity + F_contact
        % 浮力/重力: (rho_s - rho_f) * V * g
        rho_f_approx = 1000; % 假设流体密度，或从 Fluid 结构获取
        F_buoy_x = 0;
        F_buoy_y = (Particle.rho(p) - rho_f_approx) * Particle.V(p) * g_vec(2);

        Total_Fx = F_hydro_x(p) + F_buoy_x + F_cont_x;
        Total_Fy = F_hydro_y(p) + F_buoy_y + F_cont_y;
        Total_T  = T_hydro(p) + T_cont;

        % 更新速度 (Explicit Euler for sub-step)
        % 论文推荐 Crank-Nicolson，但若 dt_sub 足够小，Euler 亦可接受且实现简单
        Particle.u(p) = Particle.u(p) + (Total_Fx / m) * dt_sub;
        Particle.v(p) = Particle.v(p) + (Total_Fy / m) * dt_sub;
        Particle.omega(p) = Particle.omega(p) + (Total_T / I) * dt_sub;

        % 更新位置
        Particle.x_c(p) = Particle.x_c(p) + Particle.u(p) * dt_sub;
        Particle.y_c(p) = Particle.y_c(p) + Particle.v(p) * dt_sub;

    end % end particle loop
end % end sub-step loop

end

%% 辅助函数
function [kn, eta] = calc_coeffs(m_eff, e, T_col)
% 根据 Paper Eq. 18 计算弹簧刚度和阻尼
% e: 恢复系数
if T_col < 1e-9, T_col = 1e-9; end
ln_e = log(e);
factor = sqrt(pi^2 + ln_e^2);
kn = m_eff * (factor / T_col)^2;
eta = -2 * m_eff * ln_e / T_col;
end

function res = cross_product_2d(omega_scalar, r_vec)
% 计算 omega x r (2D平面)
% omega 沿 z 轴, r 在 xy 平面
% [0, 0, w] x [rx, ry, 0] = [-w*ry, w*rx, 0]
res = [-omega_scalar * r_vec(2); omega_scalar * r_vec(1)];
end

function lambda = lubrication_factor_pp(eps, eps_sigma, eps_dx)
% 粒子-粒子 润滑修正因子 (Eq. 32)
% eps: 归一化间隙
% 截断逻辑: if eps < eps_sigma -> use eps_sigma
%           if eps > eps_dx    -> 0 (handled outside)

val_eps = max(eps, eps_sigma);

% 理论值 (Brenner 1961 asymptotic)
lambda_raw = 1/(2*val_eps); % 主导项 O(1/eps)

% 减去 IBM 已经能解析的部分 (Correction term)
% 论文简化处理： lambda(eps) - lambda(eps_dx)
lambda_cutoff = 1/(2*eps_dx);

lambda = lambda_raw - lambda_cutoff;
if lambda < 0, lambda = 0; end
end

function lambda = lubrication_factor_pw(eps, eps_sigma, eps_dx)
% 粒子-墙 润滑修正因子 (Eq. 33)
val_eps = max(eps, eps_sigma);
lambda_raw = 1/val_eps; % 主导项 for wall is 1/eps
lambda_cutoff = 1/eps_dx;

lambda = lambda_raw - lambda_cutoff;
if lambda < 0, lambda = 0; end
end
