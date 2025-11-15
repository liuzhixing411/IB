clc;
clear;
close all;
%% Variable
L=0.1;
W=0.1;
h=0.001;
timer=0;
T=600;
M=L/h;
N=W/h;
ghostnum=3;
start=ghostnum+1;
endy=M+ghostnum;
endx=N+ghostnum;
Ny=ghostnum*2+M;
Nx=ghostnum*2+N;
X=linspace(0,W,N);
Y=linspace(0,L,M);
[XX,YY]=meshgrid(X,Y);

%% grid
grid.M=M;
grid.N=N;
grid.ghostnum=ghostnum;
grid.start=start;
grid.endy=endy;
grid.endx=endx;
grid.Ny=Ny;
grid.Nx=Nx;
grid.h=h;
%% Particle Initialization
Particle.x_c=0.05;
Particle.y_c=0.075;
Particle.r=0.005;
Particle.rho=7870;
Particle.V=pi*(Particle.r)^2; 
Particle.m=Particle.rho*Particle.V;
Particle.u=0;
Particle.v=0;
Particle.omega=0;  
Particle.I=0.5*Particle.r^2*Particle.m; 
%% Fluid初始化
Fluid.rhol=1000;
Fluid.rhog=13;
Fluid.sigma=7.28e-2;
Fluid.mul=1e-3;
Fluid.mug=1.5e-5;
% Initial Velocity
u_f_new=zeros(Ny+1,Nx+1);
v_f_new=zeros(Ny+1,Nx+1);
% Initial VOF
alpha_new=zeros(Ny,Nx);
frac=zeros(M,N);
x_c2=0.05;
y_c2=0.025;
r=0.01;
for j=1:M
    for i=1:N
        x=i*h-h/2;
        y=j*h-h/2; 
        distance=sqrt((x-x_c2)^2+(y-y_c2)^2)-r;
        if(abs(distance)< h/sqrt(2))
            count=0;
            res_y=100;
            res_x=100;
            xmin=i*h-h;
            xmax=i*h;
            ymin=j*h-h;
            ymax=j*h;
            dx=(xmax-xmin)/res_x;
            dy=(ymax-ymin)/res_y;
            for p=1:res_y
                for q=1:res_x
                    xx=xmin+q*dx-dx/2;
                    yy=ymin+p*dy-dy/2;
                    dis=sqrt((xx-x_c2)^2+(yy-y_c2)^2)-r;
                    if(dis>0)
                        count=count+1;
                    end
                end
            end

            frac(j,i)=count/(res_y*res_x);
        elseif(distance>0)
            frac(j,i)=1;

        else
            frac(j,i)=0;

        end

    end
end

alpha_new(start:endy,start:endx)=frac;
% Boundary Condition of VOF
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

% Initial pressure
% p0 
p0 = zeros(Ny,Nx); 
g_init = 9.8;
rhol_init = Fluid.rhol;
p_ref = 0.0;
for j = start:endy
    for i = start:endx
        depth = (endy - j) * h;
        p0(j,i) = p_ref + rhol_init * g_init * depth;
    end
end
% Pressure Boundary Condition (Neumann for hydrostatic)
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
% Fluid.p 存储的是动压 p' (paper eq 2)
Fluid.p = zeros(Ny,Nx); 

%% Main Loop
%% GIF
save_gif = true; 
gif_filename = 'fluid_particle_collision.gif';
figure('Position', [100, 100, 1000, 800]);

for n=1:T
    fprintf('Step %d/%d, Time = %.6f\n', n, T, timer);
    
    % Dynamic Time Step (CFL)
    eps=0.1;
    Umax=max(max(abs(Fluid.u)));
    Vmax=max(max(abs(Fluid.v)));
    Max=max(Umax,Vmax);
    if Max < 1e-10
        dt = 1e-4;  
    else
        dt=0.1*min(h/Max,eps);
    end
    timer=timer+dt;
    grid.dt=dt;
    
    % 1. 预测步 (Prediction)
    % 计算 \hat{u} (不含压力梯度和IBM力) - 对应论文 Eq (3)
    [u_pred,v_pred] = prediction(Fluid, grid);
    
    % 2. 计算IBM相互作用力
    % f_x, f_y 是加速度 (force per unit mass) - 对应论文 f
    [f_x,f_y] = Interaction(u_pred,v_pred,Particle,grid);
    
    % 3. 更新粒子 
    % 包含流体动力、浮力和重力 - 对应论文 Eq (7)
    Particle = update_particle(f_x,f_y,Particle,Fluid,grid);
    
    %  4. 更新VOF 
    alpha_new = update_alpha(Fluid,grid);
    
    % --- 5. 压力修正和流场更新 ---
    
    % 求解 Poisson 
    Fluid = update_fluid(u_pred, v_pred, alpha_new, f_x, f_y, Fluid, grid);

    %% 可视化
    
     if mod(n,5)==0  % 每5步绘制一次以提高性能
        clf;

        % 提取内部区域的速度场
        u_plot = 0.5*(Fluid.u(start:endy, start:endx)+Fluid.u(start:endy, start+1:endx+1));
        v_plot = 0.5*(Fluid.v(start:endy, start:endx)+Fluid.v(start+1:endy+1, start:endx));

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

        contour(XX,YY,Fluid.alpha(start:endy,start:endx),[0.5 0.5],'k','LineWidth',2);


        % 设置图形属性
        axis equal;
        axis([0 W 0 L]);
        title(sprintf('Fluid-Particle Interaction (Time = %.3f s, Frame = %d)', timer, n), ...
              'FontSize', 12, 'FontWeight', 'bold');
        xlabel('x (m)', 'FontSize', 11);
        ylabel('y (m)', 'FontSize', 11);

        % 添加信息文本
        info_str = {sprintf('Particle Position: (%.3f, %.3f)', Particle.x_c, Particle.y_c), ...
                   sprintf('Particle Velocity: (%.3f, %.3f) m/s', Particle.u, Particle.v), ...
                   sprintf('Max Velocity: %.4f m/s', max(speed(:)))};
        text(0.02, 0.98, info_str, 'Units', 'normalized', ...
             'VerticalAlignment', 'top', 'FontSize', 10, ...
             'BackgroundColor', 'red', 'EdgeColor', 'black');

        % 设置坐标轴方向
        set(gca, 'YDir', 'normal');

        drawnow;

        % 保存为GIF - 优化版本
        if save_gif
            % 捕获帧并调整大小来减小文件
            frame = getframe(gcf);
            im = frame2im(frame);

            
            scale_factor = 0.6; % 调整这个值来平衡质量和大小
            im_small = imresize(im, scale_factor);

            [imind, cm] = rgb2ind(im_small, 64); 

           
            if mod(n,10)==0 || n <= 10 || Particle.y_c - Particle.r <= 1.5*h
                if n == 5  % 第一次写入
                    imwrite(imind, cm, gif_filename, 'gif', ...
                           'Loopcount', inf, 'DelayTime', 0.15 );
                else
                    imwrite(imind, cm, gif_filename, 'gif', ...
                           'WriteMode', 'append', 'DelayTime', 0.15);
                end
            end
        end
    end

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

            % 应用相同的优化
            im_small = imresize(im, 0.6);
            [imind, cm] = rgb2ind(im_small, 64);
            imwrite(imind, cm, gif_filename, 'gif', ...
                   'WriteMode', 'append', 'DelayTime', 1.0, 'Compression', 'lzw');
        end
        break;
    end

end    
