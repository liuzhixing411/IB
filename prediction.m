function [u_pred, v_pred] = prediction(Fluid, grid,p0)
% Inputs from Fluid and grid structs
ghostnum = grid.ghostnum;
start = grid.start;
endy  = grid.endy;
endx  = grid.endx;
Ny = grid.Ny;
Nx = grid.Nx;
h = grid.h;
dt = grid.dt;
u_in = Fluid.u;
v_in = Fluid.v;
p_in = Fluid.p+p0;
alpha = Fluid.alpha;
rhol = Fluid.rhol; rhog = Fluid.rhog;
mul  = Fluid.mul;  mug  = Fluid.mug;
sigma = Fluid.sigma;
g = 9.8;
eps_div = 1e-12;

% build rho_cc and mu_cc (cell-centered) based on alpha
rho_cc = rhol .* alpha + rhog .* (1 - alpha);
mu_cc  = mul  .* alpha + mug  .* (1 - alpha);

%% Step1: 计算 u_pred
% 在u-face上插值得到密度和粘度
muhalfx=zeros(Ny+1,Nx+1);
rhohalfx=zeros(Ny+1,Nx+1);
for j=start:endy+1
    for i=start:endx+1
        muhalfx(j,i)=(mu_cc(j,i)+mu_cc(j,i-1))/2;
        rhohalfx(j,i)=(rho_cc(j,i)+rho_cc(j,i-1))/2;
    end
end

% 在v-face上插值得到密度和粘度
muhalfy=zeros(Ny+1,Nx+1);
rhohalfy=zeros(Ny+1,Nx+1);
for j=start:endy+1
    for i=start:endx+1
        muhalfy(j,i)=(mu_cc(j,i)+mu_cc(j-1,i))/2;
        rhohalfy(j,i)=(rho_cc(j,i)+rho_cc(j-1,i))/2;
    end
end

% 对流项
adux=zeros(Ny+1,Nx+1);
aduy=zeros(Ny+1,Nx+1);

for j=start:endy+1
    for i=start:endx+1
        if(u_in(j,i)>0)
            adux(j,i)=u_in(j,i)*(3*u_in(j,i)-4*u_in(j,i-1)+u_in(j,i-2))/(2*h);
        else
            adux(j,i)=u_in(j,i)*(-3*u_in(j,i)+4*u_in(j,i+1)-u_in(j,i+2))/(2*h);
        end
        
        v_at_u_node = 0.25 * (v_in(j,i-1) + v_in(j,i) + v_in(j+1,i-1) + v_in(j+1,i));
        if(v_at_u_node>0)
            aduy(j,i)=v_at_u_node*(3*u_in(j,i)-4*u_in(j-1,i)+u_in(j-2,i))/(2*h);
        else
            aduy(j,i)=v_at_u_node*(-3*u_in(j,i)+4*u_in(j+1,i)-u_in(j+2,i))/(2*h);
        end
    end
end
adu=adux+aduy;

% 扩散项的速度梯度
diffux=zeros(Ny+1,Nx+1);
diffuy=zeros(Ny+1,Nx+1);
diffvx=zeros(Ny+1,Nx+1);
diffvy=zeros(Ny+1,Nx+1);

for j=start:endy+1
    for i=start:endx+1
        if(u_in(j,i)>0)
            diffux(j,i)=(u_in(j,i)-u_in(j,i-1))/h;
            diffvx(j,i)=(v_in(j,i)-v_in(j,i-1))/h;
        else
            diffux(j,i)=(u_in(j,i+1)-u_in(j,i))/h;
            diffvx(j,i)=(v_in(j,i+1)-v_in(j,i))/h;
        end

        if(v_in(j,i)>0)
            diffuy(j,i)=(u_in(j,i)-u_in(j-1,i))/h;
            diffvy(j,i)=(v_in(j,i)-v_in(j-1,i))/h;
        else
            diffuy(j,i)=(u_in(j+1,i)-u_in(j,i))/h;
            diffvy(j,i)=(v_in(j+1,i)-v_in(j,i))/h;
        end
    end  
end

% 计算u方向的扩散项
diffu=zeros(Ny+1,Nx+1);
for j=start:endy+1
    for i=start:endx+1
        % Term 1: div(mu * grad(u)) - Laplacian part
        laplacian_u = (u_in(j+1,i) - 2*u_in(j,i) + u_in(j-1,i))/h^2 + ...
                      (u_in(j,i+1) - 2*u_in(j,i) + u_in(j,i-1))/h^2;
        diff1 = laplacian_u * muhalfx(j,i);
        
        % Term 2: d(mu)/dx * du/dx
        dmu_dx = (mu_cc(j,i) - mu_cc(j,i-1))/h;
        diff2x = dmu_dx * 2 * diffux(j,i);
        
        % Term 3: d(mu)/dy * (du/dy + dv/dx)
        dmu_dy_N = (mu_cc(j+1,i) + mu_cc(j+1,i-1)) * 0.5;
        dmu_dy_S = (mu_cc(j-1,i) + mu_cc(j-1,i-1)) * 0.5;
        dmu_dy = (dmu_dy_N - dmu_dy_S) / (2*h);
        
        % 更准确的dv/dx在u-face处
        dvdx_at_uface = (diffvx(j,i) + diffvx(j+1,i))/2;
        diff2y = dmu_dy * (diffuy(j,i) + dvdx_at_uface);
        
        diffu(j,i) = (diff1 + diff2x + diff2y) / (rhohalfx(j,i) + eps_div);
    end
end

% 表面张力 (CSF model)
alpha_x = zeros(Ny, Nx); alpha_y = zeros(Ny, Nx);
for j = 2:Ny-1
    for i = 2:Nx-1
        alpha_x(j,i) = (alpha(j, i+1) - alpha(j, i-1)) / (2*h);
        alpha_y(j,i) = (alpha(j+1, i) - alpha(j-1, i)) / (2*h);
    end
end

% 曲率计算
kappa = zeros(Ny, Nx);
for j = 2:Ny-1
    for i = 2:Nx-1
        grad_abs = sqrt(alpha_x(j,i)^2 + alpha_y(j,i)^2) + eps_div;
        
        % 计算法向量散度 - 修正版
        % 在每个点重新计算grad_abs
        grad_abs_ip1 = sqrt(alpha_x(j,i+1)^2 + alpha_y(j,i+1)^2) + eps_div;
        grad_abs_im1 = sqrt(alpha_x(j,i-1)^2 + alpha_y(j,i-1)^2) + eps_div;
        grad_abs_jp1 = sqrt(alpha_x(j+1,i)^2 + alpha_y(j+1,i)^2) + eps_div;
        grad_abs_jm1 = sqrt(alpha_x(j-1,i)^2 + alpha_y(j-1,i)^2) + eps_div;
        
        div_n = (alpha_x(j,i+1)/grad_abs_ip1 - alpha_x(j,i-1)/grad_abs_im1)/(2*h) + ...
                (alpha_y(j+1,i)/grad_abs_jp1 - alpha_y(j-1,i)/grad_abs_jm1)/(2*h);
        kappa(j,i) = -div_n;
    end
end

% 边界外推
for j=1:ghostnum
    for i=1:Nx
        kappa(j,i)=kappa(2*ghostnum+1-j,i);
    end
end
for j=endy+1:Ny
    for i=1:Nx
        kappa(j,i)=kappa(2*endy+1-j,i);
    end
end
for j=1:Ny
    for i=1:ghostnum
        kappa(j,i)=kappa(j,2*ghostnum+1-i);
    end
end
for j=1:Ny
    for i=endx+1:Nx
        kappa(j,i)=kappa(j,2*endx+1-i);
    end
end

% CSF force at cell centers
fcell_x = sigma .* kappa .* alpha_x;
fcell_y = sigma .* kappa .* alpha_y;

% 插值到面上
surfu = zeros(Ny+1, Nx+1);
surfv = zeros(Ny+1, Nx+1);
for jf = 1:Ny+1
    for ifc = 1:Nx+1
        jc = max(min(jf,Ny),1);
        iL = max(min(ifc-1,Nx),1);
        iR = max(min(ifc,Nx),1);
        surfu(jf,ifc) = 0.5 * (fcell_x(jc, iL) + fcell_x(jc, iR)) / (rhohalfx(jf,ifc) + eps_div);
        
        jrB = max(min(jf-1,Ny),1);
        jrT = max(min(jf,Ny),1);
        ic = max(min(ifc,Nx),1);
        surfv(jf,ifc) = 0.5 * (fcell_y(jrB, ic) + fcell_y(jrT, ic)) / (rhohalfy(jf,ifc) + eps_div);
    end
end

% **修正：压力梯度计算 - 分别处理u和v面**
Pu = zeros(Ny+1, Nx+1);
Pv = zeros(Ny+1, Nx+1);

% u-face上的压力梯度 (从cell center到u-face)
for j = start:endy+1
    for i = start:endx+1
        dpdx = (p_in(j, i) - p_in(j, i-1)) / h;
        Pu(j,i) = dpdx / (rhohalfx(j,i) + eps_div);
    end
end

% v-face上的压力梯度 (从cell center到v-face)
for j = start:endy+1
    for i = start:endx+1
        dpdy = (p_in(j, i) - p_in(j-1, i)) / h;
        Pv(j,i) = dpdy / (rhohalfy(j,i) + eps_div);
    end
end

% 求解u_pred (包含压力梯度项)
u_pred=u_in+(-adu+diffu+surfu - Pu)*dt; 

% 边界条件
for j=1:Ny+1
    for i=1:ghostnum
        u_pred(j,i)=-u_pred(j,2*ghostnum+2-i);
    end
end
for j=1:Ny+1
    for i=endx+1:Nx+1
        u_pred(j,i)=-u_pred(j,2*endx+2-i);
    end
end
for j=1:ghostnum
    for i=1:Nx+1
        u_pred(j,i)=-u_pred(2*ghostnum+1-j,i);
    end
end
for j=endy+1:Ny+1
    for i=1:Nx+1
        u_pred(j,i)=-u_pred(2*endy+1-j,i);
    end
end
u_pred(:,ghostnum+1)=0.0;
u_pred(:,endx+1)=0.0;

%% Step2: 计算v_pred
% 对流项
advx=zeros(Ny+1,Nx+1);
advy=zeros(Ny+1,Nx+1);

for j=start:endy+1
    for i=start:endx+1
        u_at_v_node = 0.25 * (u_in(j-1,i) + u_in(j-1,i+1) + u_in(j,i) + u_in(j,i+1));
        if(u_at_v_node>0)
            advx(j,i)=u_at_v_node*(3*v_in(j,i)-4*v_in(j,i-1)+v_in(j,i-2))/(2*h);
        else
            advx(j,i)=u_at_v_node*(-3*v_in(j,i)+4*v_in(j,i+1)-v_in(j,i+2))/(2*h);
        end
        
        if(v_in(j,i)>0)
            advy(j,i)=v_in(j,i)*(3*v_in(j,i)-4*v_in(j-1,i)+v_in(j-2,i))/(2*h);
        else
            advy(j,i)=v_in(j,i)*(-3*v_in(j,i)+4*v_in(j+1,i)-v_in(j+2,i))/(2*h);
        end
    end
end
adv=advx+advy;

% 计算v方向的扩散项
diffv=zeros(Ny+1,Nx+1);
for j=start:endy+1
    for i=start:endx+1
        % Term 1: div(mu * grad(v))
        laplacian_v = (v_in(j+1,i) - 2*v_in(j,i) + v_in(j-1,i))/h^2 + ...
                      (v_in(j,i+1) - 2*v_in(j,i) + v_in(j,i-1))/h^2;
        diff1 = muhalfy(j,i) * laplacian_v;
        
        % Term 2: d(mu)/dy * dv/dy
        dmu_dy = (mu_cc(j,i) - mu_cc(j-1,i))/h;
        diff2y = dmu_dy * 2 * diffvy(j,i);
        
        % Term 3: d(mu)/dx * (du/dy + dv/dx)
        dmu_dx_E = (mu_cc(j,i+1) + mu_cc(j-1,i+1)) * 0.5;
        dmu_dx_W = (mu_cc(j,i-1) + mu_cc(j-1,i-1)) * 0.5;
        dmu_dx = (dmu_dx_E - dmu_dx_W) / (2*h);
        
        % 更准确的du/dy在v-face处
        dudy_at_vface = (diffuy(j,i) + diffuy(j,i+1))/2;
        diff2x = dmu_dx * (dudy_at_vface + diffvx(j,i));
        
        diffv(j,i) = (diff1 + diff2x + diff2y) / (rhohalfy(j,i) + eps_div);
    end
end

% 求解v_pred (包含压力梯度和重力)
v_pred=v_in+(-adv+diffv+surfv - Pv - g)*dt;
% 边界条件
for j=1:Ny+1
    for i=1:ghostnum
        v_pred(j,i)=-v_pred(j,2*ghostnum+1-i);
    end
end
for j=1:Ny+1
    for i=endx+1:Nx+1
        v_pred(j,i)=-v_pred(j,2*endx+1-i);
    end
end
for j=1:ghostnum
    for i=1:Nx+1
        v_pred(j,i)=-v_pred(2*ghostnum+2-j,i);
    end
end
for j=endy+1:Ny+1
    for i=1:Nx+1
        v_pred(j,i)=-v_pred(2*endy+2-j,i);
    end
end
v_pred(ghostnum+1,:)=0.0;
v_pred(endy+1,:)=0.0;

end