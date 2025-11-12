function pcc = Poisson(rho_cc, u_new, v_new, u_old, v_old, dx, dy, dt, ghostnum, f_x, f_y)

[Ny, Nx] = size(rho_cc);
start = ghostnum + 1;
endy  = Ny - ghostnum;
endx  = Nx - ghostnum;
M = endy - start + 1;
N = endx - start + 1;
nCells = M * N;

% beta at cell centers
beta_cc = 1 ./ rho_cc;

% allocate assembly arrays
I = zeros(5*nCells,1); J = I; S = I; ptr = 0;
b = zeros(nCells,1);

idx = @(jj,ii) (ii - start) * M + (jj - start) + 1;

small = 1e-14;

% -------------------- 对流项计算 --------------------
adux = zeros(Ny+1,Nx+1); aduy = zeros(Ny+1,Nx+1);
advx = zeros(Ny+1,Nx+1); advy = zeros(Ny+1,Nx+1);

% u方向的对流加速度
for j=start:endy+1
    for i=start:endx+1
        if(u_new(j,i)>0)
            adux(j,i)=u_new(j,i)*(3*u_new(j,i)-4*u_new(j,i-1)+u_new(j,i-2))/(2*dx);
        else
            adux(j,i)=u_new(j,i)*(-3*u_new(j,i)+4*u_new(j,i+1)-u_new(j,i+2))/(2*dx);
        end
        
        v_at_u_node = 0.25 * (v_new(j,i-1) + v_new(j,i) + v_new(j+1,i-1) + v_new(j+1,i));
        if(v_at_u_node>0)
            aduy(j,i)=v_at_u_node*(3*u_new(j,i)-4*u_new(j-1,i)+u_new(j-2,i))/(2*dy);
        else
            aduy(j,i)=v_at_u_node*(-3*u_new(j,i)+4*u_new(j+1,i)-u_new(j+2,i))/(2*dy);
        end
    end
end
adu=adux+aduy;

% v方向的对流加速度
for j=start:endy+1
    for i=start:endx+1
        u_at_v_node = 0.25 * (u_new(j-1,i) + u_new(j-1,i+1) + u_new(j,i) + u_new(j,i+1));
        % **修正：v的x方向导数应该除以dx**
        if(u_at_v_node>0)
            advx(j,i)=u_at_v_node*(3*v_new(j,i)-4*v_new(j,i-1)+v_new(j,i-2))/(2*dx);
        else
            advx(j,i)=u_at_v_node*(-3*v_new(j,i)+4*v_new(j,i+1)-v_new(j,i+2))/(2*dx);
        end
        
        if(v_new(j,i)>0)
            advy(j,i)=v_new(j,i)*(3*v_new(j,i)-4*v_new(j-1,i)+v_new(j-2,i))/(2*dy);
        else
            advy(j,i)=v_new(j,i)*(-3*v_new(j,i)+4*v_new(j+1,i)-v_new(j+2,i))/(2*dy);
        end
    end
end
adv=advx+advy;

% -------------------- 对流加速度散度 div((u·∇)u) --------------------
div_a_new = zeros(Ny,Nx);

for j = start:endy
    for i = start:endx
        % u方向对流加速度的x分量散度
        dax_dx_new = (adu(j, i+1) - adu(j, i)) / dx;
        % v方向对流加速度的y分量散度
        day_dy_new = (adv(j+1, i) - adv(j, i)) / dy;
        div_a_new(j,i) = dax_dx_new + day_dy_new;
    end
end

% -------------------- 外力散度 div(f) --------------------
div_f = zeros(Ny,Nx);
for j = start:endy
    for i = start:endx
        div_f(j,i) = (f_x(j, i+1) - f_x(j, i)) / dx + (f_y(j+1, i) - f_y(j, i)) / dy;
    end
end

% -------------------- 速度散度 (不可压缩性残差) --------------------
divU_new = zeros(Ny,Nx);
divU_old = zeros(Ny,Nx);
for j = start:endy
    for i = start:endx
        divU_new(j,i) = (u_new(j, i+1) - u_new(j, i)) / dx + (v_new(j+1, i) - v_new(j, i)) / dy;
        divU_old(j,i) = (u_old(j, i+1) - u_old(j, i)) / dx + (v_old(j+1, i) - v_old(j, i)) / dy;
    end
end

% 时间导数项 dD/dt ≈ (divU_new - divU_old)/dt
dD_dt = (divU_new - divU_old) / dt;

% ---------------------------------------------------------------------
% 组装矩阵 A (变系数) 和右端项 b = -div_a_new - dD_dt + div_f
% ---------------------------------------------------------------------
for j = start:endy
    for i = start:endx
        k = idx(j,i);

        % 面上的调和平均beta (1/rho)
        beta_e = 2 * beta_cc(j,i)   * beta_cc(j,i+1) / (beta_cc(j,i)   + beta_cc(j,i+1)   + small);
        beta_w = 2 * beta_cc(j,i-1) * beta_cc(j,i)   / (beta_cc(j,i-1) + beta_cc(j,i)     + small);
        beta_n = 2 * beta_cc(j+1,i) * beta_cc(j,i)   / (beta_cc(j+1,i) + beta_cc(j,i)     + small);
        beta_s = 2 * beta_cc(j,i)   * beta_cc(j-1,i) / (beta_cc(j,i)   + beta_cc(j-1,i)   + small);

        aE = beta_e / dx^2;
        aW = beta_w / dx^2;
        aN = beta_n / dy^2;
        aS = beta_s / dy^2;

        diagA = (aE + aW + aN + aS);

        % 右端项 S = -div_a - dD_dt + div_f
        Sval = - div_a_new(j,i) - dD_dt(j,i) + div_f(j,i);
        b(k) = -Sval;

        % West neighbor
        if i-1 >= start
            ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j,i-1); S(ptr)=-aW;
        else
            diagA = diagA - aW;
        end

        % East neighbor
        if i+1 <= endx
            ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j,i+1); S(ptr)=-aE;
        else
            diagA = diagA - aE;
        end

        % South neighbor
        if j-1 >= start
            ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j-1,i); S(ptr)=-aS;
        else
            diagA = diagA - aS;
        end

        % North neighbor
        if j+1 <= endy
            ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j+1,i); S(ptr)=-aN;
        else
            diagA = diagA - aN;
        end

        % diagonal
        ptr = ptr + 1; I(ptr)=k; J(ptr)=k; S(ptr)=diagA;
    end
end

% 形成稀疏矩阵
I = I(1:ptr); J = J(1:ptr); S = S(1:ptr);
A = sparse(I,J,S,nCells,nCells);

% 固定参考压力 (中心点)
ref_i = round((start+endx)/2);
ref_j = round((start+endy)/2);
ref_k = idx(ref_j, ref_i);
A(ref_k, :) = 0;
A(ref_k, ref_k) = 1;
b(ref_k) = 0;

% 求解
pvec = A \ b;

% 重塑为带ghost的pcc
pcc = zeros(Ny, Nx);
p_interior = reshape(pvec, M, N);
pcc(start:endy, start:endx) = p_interior;

% 镜像边界条件
for j = 1:ghostnum
    mirror_j = 2*ghostnum + 1 - j;
    pcc(j, :) = pcc(mirror_j, :);
end
for j = endy+1:Ny
    mirror_j = 2*endy + 1 - j;
    pcc(j, :) = pcc(mirror_j, :);
end
for i = 1:ghostnum
    mirror_i = 2*ghostnum + 1 - i;
    pcc(:, i) = pcc(:, mirror_i);
end
for i = endx+1:Nx
    mirror_i = 2*endx + 1 - i;
    pcc(:, i) = pcc(:, mirror_i);
end

end