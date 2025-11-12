clc;
clear;
close all;

%% 变量定义
% 求解区域
L=0.1;
W=0.1;
h=0.001;
timer=0;
T=1000;
M=L/h;
N=W/h;
ghostnum=3;
start=ghostnum+1;
endy=M+ghostnum;
endx=N+ghostnum;
Ny=ghostnum*2+M;
Nx=ghostnum*2+N;

% 网格
X=linspace(0,W,N);
Y=linspace(0,L,M);
[X,Y]=meshgrid(X,Y);

% 尺寸结构体
grid.M=M;
grid.N=N;
grid.ghostnum=ghostnum;
grid.start=start;
grid.endy=endy;
grid.endx=endx;
grid.Ny=Ny;
grid.Nx=Nx;
grid.h=h;

%% Particle初始化
Particle.x_c=0.05;
Particle.y_c=0.075;
Particle.r=0.01;
Particle.rho=7870;
Particle.V=pi*(Particle.r)^2;
Particle.m=Particle.rho*Particle.V;
Particle.u=0;
Particle.v=0;
Particle.omega=0;  % **修正：omega应该是标量，不是数组**
Particle.I=0.5*Particle.r^2*Particle.m;


%% Fluid初始化
Fluid.rhol=1000;
Fluid.rhog=1.3;
Fluid.sigma=7.28e-2;
Fluid.mul=1e-3;
Fluid.mug=1.5e-5;

% 初始速度场
u_f_new=zeros(Ny+1,Nx+1);
v_f_new=zeros(Ny+1,Nx+1);

% 初始alpha场 (全部是液体)
alpha_new=ones(Ny,Nx);

% 应用边界条件到alpha
for j=1:ghostnum
    for i=1:Nx
        alpha_new(j,i)=alpha_new(2*ghostnum+1-j,i);
    end
end
for j=endy+1:Ny
    for i=1:Nx
        alpha_new(j,i)=alpha_new(2*endy+1-j,i);
    end
end
for j=1:Ny
    for i=1:ghostnum
        alpha_new(j,i)=alpha_new(j,2*ghostnum+1-i);
    end
end
for j=1:Ny
    for i=endx+1:Nx
        alpha_new(j,i)=alpha_new(j,2*endx+1-i);
    end
end

Fluid.alpha=alpha_new;
Fluid.u=u_f_new;
Fluid.v=v_f_new;

% 初始压力场（静水压分布）
p0 = zeros(Ny,Nx);
p_ref = 0.0;
g_init = 9.8;
rhol_init = Fluid.rhol;

for j = start:endy
    for i = start:endx
        % 压力随深度线性增加: p = p_ref + ρgh
        % y = 0 在底部，所以 h = (endy - j) * grid.h
        depth = (endy - j) * h;
        p0(j,i) = p_ref + rhol_init * g_init * depth;
    end
end

% 压力边界条件
for j=1:ghostnum
    for i=1:Nx
        p0(j,i)=p0(2*ghostnum+1-j,i);
    end
end
for j=endy+1:Ny
    for i=1:Nx
        p0(j,i)=p0(2*endy+1-j,i);
    end
end
for j=1:Ny
    for i=1:ghostnum
        p0(j,i)=p0(j,2*ghostnum+1-i);
    end
end
for j=1:Ny
    for i=endx+1:Nx
        p0(j,i)=p0(j,2*endx+1-i);
    end
end

Fluid.p = zeros(Ny,Nx);

%% 主循环
%% GIF输出设置
save_gif = true; % 设置为true输出GIF，false则不输出
gif_filename = 'fluid_particle_simulation.gif';
if save_gif
    fprintf('GIF输出已启用，文件将保存为: %s\n', gif_filename);
end

%% 主循环
figure('Position', [100, 100, 1000, 800]);

for n=1:T
    fprintf('Step %d/%d, Time = %.6f\n', n, T, timer);
    
    % 动态时间步长 (CFL条件)
    eps=0.01;
    Umax=max(max(abs(Fluid.u)));
    Vmax=max(max(abs(Fluid.v)));
    Max=max(Umax,Vmax);
    if Max < 1e-10
        dt = 1e-4;  % 初始时刻的小时间步
    else
        dt=0.1*min(h/Max,eps);
    end
    timer=timer+dt;
    grid.dt=dt;
    
    [u_pred,v_pred]=prediction(Fluid,grid,p0);
    
    [f_x,f_y]=Interaction(u_pred,v_pred,Particle,grid);
    

    Particle=update_particle(f_x,f_y,Particle,grid);
    

    alpha_new=update_alpha(Fluid,grid);
    

    Fluid=update_fluid(u_pred,v_pred,alpha_new,f_x,f_y,Fluid,grid);
    
    %% 可视化
    % if mod(n,10)==0
    %     clf;
    %     hold on;
    %     axis equal;
    % 
    %     % 绘制速度场
    %     u_plot = Fluid.u(start:endy,start:endx);
    %     v_plot = Fluid.v(start:endy,start:endx);
    %     % 插值到cell center
    %     u_cc = 0.5*(u_plot + u_plot);
    %     v_cc = 0.5*(v_plot + v_plot);
    %     %contour(X,Y,sqrt(v_cc.^2+u_cc.^2),10,'b');
    % 
    %     streamslice(X,Y,u_cc,v_cc,5);
    % 
    %     % 绘制粒子
    %     rectangle('Position', [Particle.x_c-Particle.r, Particle.y_c-Particle.r, ...
    %                            2*Particle.r, 2*Particle.r], ...
    %               'Curvature', [1, 1], ...
    %               'FaceColor', 'red', ...
    %               'EdgeColor', 'red', ...
    %               'LineWidth', 2);
    % 
    %     xlim([0 W]);
    %     ylim([0 L]);
    %     title(sprintf('Time = %.4f s, Step = %d', timer, n));
    %     xlabel('x (m)');
    %     ylabel('y (m)');
    %     drawnow;
    % end
    % 
    % % 检查粒子是否到达底部
    % if Particle.y_c - Particle.r <= 0
    %     fprintf('Particle reached bottom at step %d\n', n);
    %     break;
    % end
    %% 新的可视化 - 彩色速度分布
    if mod(n,5)==0  % 每5步绘制一次以提高性能
        clf;
        
        % 提取内部区域的速度场
        u_plot = Fluid.u(start:endy, start:endx);
        v_plot = Fluid.v(start:endy, start:endx);
        
        % 计算速度大小
        speed = sqrt(u_plot.^2 + v_plot.^2);
        
        % 创建子图布局
        ax = axes;
        
        % 绘制彩色速度背景
        imagesc([0 W], [0 L], speed);
        hold on;
        
        % 设置颜色映射
        colormap(jet);
        colorbar;
        clim([0 max(0.01, max(speed(:)))]); % 动态调整颜色范围
        
        % 绘制流线
        [x_grid, y_grid] = meshgrid(linspace(0, W, size(u_plot,2)), ...
                                   linspace(0, L, size(u_plot,1)));
        h_stream = streamslice(x_grid, y_grid, u_plot, v_plot, 2);
        set(h_stream, 'Color', 'white', 'LineWidth', 1);
        
        % 绘制粒子
        rectangle('Position', [Particle.x_c-Particle.r, Particle.y_c-Particle.r, ...
                               2*Particle.r, 2*Particle.r], ...
                  'Curvature', [1, 1], ...
                  'FaceColor', 'red', ...
                  'EdgeColor', 'red', ...
                  'LineWidth', 2);
        
        % 设置图形属性
        axis equal;
        axis([0 W 0 L]);
        title(sprintf('流体-粒子相互作用 (时间 = %.3f s, 步数 = %d)', timer, n), ...
              'FontSize', 12, 'FontWeight', 'bold');
        xlabel('x (m)', 'FontSize', 11);
        ylabel('y (m)', 'FontSize', 11);
        
        % 添加信息文本
        info_str = {sprintf('粒子位置: (%.3f, %.3f)', Particle.x_c, Particle.y_c), ...
                   sprintf('粒子速度: (%.3f, %.3f) m/s', Particle.u, Particle.v), ...
                   sprintf('最大流速: %.4f m/s', max(speed(:)))};
        text(0.02, 0.98, info_str, 'Units', 'normalized', ...
             'VerticalAlignment', 'top', 'FontSize', 10, ...
             'BackgroundColor', 'white', 'EdgeColor', 'black');
        
        % 设置坐标轴方向
        set(gca, 'YDir', 'normal');
        
        drawnow;
        
        % 保存为GIF
        if save_gif
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            
            if n == 5  % 第一次写入
                imwrite(imind, cm, gif_filename, 'gif', ...
                       'Loopcount', inf, 'DelayTime', 0.1);
            else
                imwrite(imind, cm, gif_filename, 'gif', ...
                       'WriteMode', 'append', 'DelayTime', 0.1);
            end
        end
    end
    
    % 检查粒子是否到达底部
    if Particle.y_c - Particle.r <= 0
        fprintf('粒子已到达底部，步数: %d\n', n);
        
        % 最终帧
        if save_gif
            % 添加最终状态的文本
            text(0.5, 0.5, '模拟完成 - 粒子已到达底部', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                 'FontSize', 14, 'FontWeight', 'bold', ...
                 'BackgroundColor', 'yellow', 'EdgeColor', 'red');
            drawnow;
            
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            imwrite(imind, cm, gif_filename, 'gif', ...
                   'WriteMode', 'append', 'DelayTime', 1.0);
        end
        break;
    end
end