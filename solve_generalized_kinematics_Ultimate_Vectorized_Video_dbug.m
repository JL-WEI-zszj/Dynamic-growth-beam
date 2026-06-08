function solve_generalized_kinematics_Ultimate_Vectorized_Video_dbug()

% ==================== 1. 用户自定义基础设置 ====================
h0 = 1/100;
eps = 1/1000;
L = 1;
nGrid = 500;          % 偶数高密度网格
Tmax = 20;

PlotMode = 'Relative'; 
Rayleigh_beta = 5e-8; % 物理降噪阻尼

% 【视频生成开关】
MakeVideo = true;               % 设为 true 则生成视频

currentTime = datestr(now, 'HHMMSS'); 
VideoName = sprintf('Rod_Evolution_%s.mp4', currentTime);

% ==================== 2. 自定义基座运动与预拉伸场 ====================
syms t_sym S_sym real
fprintf('========== 开始四阶高精度极速仿真 (0406边界修正版 + 视频生成) ==========\n');
fprintf('正在编译符号运算函数...\n');

cx_sym = 0 * S_sym + 0 * t_sym; 
cz_sym = 0 * S_sym + 0 * t_sym; 

Nl10_sym = 1 + 0 * t_sym;   
Nl11_sym = sin(t_sym);  % 用户修改     
Nl12_sym = 0 * S_sym + 0 * t_sym;

cxddot_sym = diff(cx_sym, t_sym, 2); czddot_sym = diff(cz_sym, t_sym, 2);
Nl10_S_sym = diff(Nl10_sym, S_sym); Nl10_t_sym = diff(Nl10_sym, t_sym); Nl10_tt_sym = diff(Nl10_t_sym, t_sym);
Nl11_S_sym = diff(Nl11_sym, S_sym); Nl11_t_sym = diff(Nl11_sym, t_sym);

F = struct();
F.cx = matlabFunction(cx_sym, 'Vars', {t_sym}); F.cz = matlabFunction(cz_sym, 'Vars', {t_sym});
F.cxddot = matlabFunction(cxddot_sym, 'Vars', {t_sym}); F.czddot = matlabFunction(czddot_sym, 'Vars', {t_sym});
F.Nl10 = matlabFunction(Nl10_sym, 'Vars', {S_sym, t_sym}); F.Nl10_S = matlabFunction(Nl10_S_sym, 'Vars', {S_sym, t_sym});
F.Nl10_t = matlabFunction(Nl10_t_sym, 'Vars', {S_sym, t_sym}); F.Nl10_tt = matlabFunction(Nl10_tt_sym,'Vars', {S_sym, t_sym});
F.Nl11 = matlabFunction(Nl11_sym, 'Vars', {S_sym, t_sym}); F.Nl11_S = matlabFunction(Nl11_S_sym, 'Vars', {S_sym, t_sym});
F.Nl11_t = matlabFunction(Nl11_t_sym, 'Vars', {S_sym, t_sym});
F.Nl12 = matlabFunction(Nl12_sym, 'Vars', {S_sym, t_sym});

% ==================== 3. 核心求解过程 ====================
fprintf('编译完成。开始基于正定质量矩阵的向量化高速 PDE 求解...\n');
tic; 
[t, Y, S] = solvePDE_MassMatrix_Vectorized(h0, eps, nGrid, Tmax, F, Rayleigh_beta);
nNodes = length(S); g_sol = Y(:, 1:nNodes); g_t = Y(:, nNodes+1:end);  
time_pde = toc;

fprintf('PDE 求解完成。共步进 %d 个时间点，耗时 %.2f 秒。\n', length(t), time_pde);
fprintf('正在计算动态张量及坐标映射...\n');

[~, ~, f_vals, x0_vals, z0_vals] = ...
    compute_kinematics_PostProcess(t, S, g_sol, g_t, h0, eps, F, PlotMode);

fprintf('计算完毕，开始生成图表和动画...\n');

% ==================== 4. 原有静态结果可视化 ====================

% 图1：自由端 Z 坐标 (保持不变)
figure('Name', '自由端 Z 坐标', 'Position', [600, 100, 500, 400]);
plot(t, z0_vals(:, end), 'LineWidth', 2, 'Color', '#D95319');
xlabel('时间 \tau'); ylabel('Z_0(1,\tau) (m)'); title('自由端 Z 向位移演化'); grid on;

% ==================== 4.5. 输出末端 Z 位移 CSV 文件 ====================
z_tip = z0_vals(:, end);                          
output_data = [t(:), z_tip(:)];                   
csv_filename = sprintf('Rod_TipZ_Displacement_%s.csv', currentTime);
fid = fopen(csv_filename, 'w');
fprintf(fid, 'Time (s),Tip_Z_Displacement (m)\n');
fclose(fid);
writematrix(output_data, csv_filename, 'WriteMode', 'append');
fprintf('末端 Z 位移数据已保存至: %s\n', csv_filename);

time_slices = 0:2:Tmax; 
colors = turbo(length(time_slices)); 

% 图2：空间形变姿态叠加
figure('Name', '空间形变姿态叠加', 'Position', [600, 550, 500, 400]); hold on; 

box off; 

for k = 1:length(time_slices)
    [~, idx] = min(abs(t - time_slices(k)));
    plot(x0_vals(idx,:), z0_vals(idx,:), 'LineWidth', 2, 'Color', colors(k,:), 'DisplayName', sprintf('t = %d s', time_slices(k)));
    plot(x0_vals(idx,1), z0_vals(idx,1), 'o', 'Color', colors(k,:), 'MarkerSize', 6, 'HandleVisibility','off');
    plot(x0_vals(idx,end), z0_vals(idx,end), '^', 'Color', colors(k,:), 'MarkerSize', 6, 'MarkerFaceColor', colors(k,:), 'HandleVisibility','off');
end

title(''); 
xlabel('$\overline{x}^{(0)}$ (m)', 'Interpreter', 'latex', 'FontSize', 12); 
ylabel('$\overline{z}^{(0)}$ (m)', 'Interpreter', 'latex', 'FontSize', 12); 
legend('Location', 'best'); 
grid off; 
ax = gca;
x_min = min(x0_vals(:));
x_max = max(x0_vals(:));

if x_min >= 0
    n_neg = 0; n_pos = 5;
elseif x_max <= 0
    n_neg = 5; n_pos = 0;
else
    best_dx = inf; best_n = 1;
    for n = 1:4
        temp_dx = max(abs(x_min)/n, x_max/(5-n));
        if temp_dx < best_dx
            best_dx = temp_dx; best_n = n;
        end
    end
    n_neg = best_n; n_pos = 5 - best_n;
end

if n_neg == 0
    min_dx = x_max / 5;
elseif n_pos == 0
    min_dx = abs(x_min) / 5;
else
    min_dx = max(abs(x_min)/n_neg, x_max/n_pos); 
end

dx = ceil(min_dx / 0.5) * 0.5;
if dx == 0; dx = 0.5; end

x_major = -n_neg*dx : dx : n_pos*dx;
x_limits = [x_major(1), x_major(end)];

xlim(x_limits);
ax.XTick = x_major;
ax.XAxis.MinorTickValues = x_major(1:end-1) + dx/2; % 副刻度放在主刻度正中间
ax.XMinorTick = 'on';

z_min = min(z0_vals(:));
z_max = max(z0_vals(:));
z_abs_max = max(abs(z_min), abs(z_max));

d = ceil(z_abs_max / 2 / 0.05) * 0.05; 
if d == 0; d = 0.05; end 

y_major = [-2*d, -d, 0, d, 2*d];
y_limits = [y_major(1), y_major(end)];

ylim(y_limits); 

ax.YTick = y_major;
ax.YAxis.MinorTickValues = y_major(1:end-1) + d/2;
ax.YMinorTick = 'on';

ax.TickDir = 'in';

plot(x_limits, [y_limits(2), y_limits(2)], '-', 'Color', ax.XColor, 'LineWidth', ax.LineWidth, 'HandleVisibility', 'off'); % 顶部线
plot([x_limits(2), x_limits(2)], y_limits, '-', 'Color', ax.YColor, 'LineWidth', ax.LineWidth, 'HandleVisibility', 'off'); % 右侧线

% ==================== 5. 动态构型演化视频生成 ====================
fps = 30;                     % 视频帧率
play_speed = 1.0;             % 播放倍速（1.0 = 真实时间）
dt_frame = play_speed / fps;
t_uniform = 0 : dt_frame : Tmax;

frame_indices = zeros(1, length(t_uniform));
for k = 1:length(t_uniform)
    [~, idx] = min(abs(t - t_uniform(k)));
    frame_indices(k) = idx;
end
frame_indices = unique(frame_indices);   % 去重

figure('Name', '杆件构型演化视频动画', 'Position', [300, 300, 700, 600]); 
hold on; box on; axis equal; grid on;

x_min = min(x0_vals(:)); x_max = max(x0_vals(:));
z_min = min(z0_vals(:)); z_max = max(z0_vals(:));
margin_x = max(0.1, (x_max - x_min)*0.1); margin_z = max(0.1, (z_max - z_min)*0.1);
xlim([x_min - margin_x, x_max + margin_x]); ylim([z_min - margin_z, z_max + margin_z]);

if strcmp(PlotMode, 'Relative')
    xlabel('相对坐标 x_0 (m)'); ylabel('相对坐标 z_0 (m)'); title('杆件随动构型动态演化 (真实时间同步)');
else
    xlabel('全局坐标 X_0 (m)'); ylabel('全局坐标 Z_0 (m)'); title('杆件全局运动轨迹动态演化 (真实时间同步)');
end

h_rod  = plot(NaN, NaN, 'LineWidth', 2.5, 'Color', '#0072BD');
h_root = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'Color', '#D95319', 'MarkerFaceColor', '#D95319');
h_tip  = plot(NaN, NaN, '^', 'MarkerSize', 8, 'Color', '#EDB120', 'MarkerFaceColor', '#EDB120');
h_time = text(x_min, z_max, '', 'FontSize', 12, 'FontWeight', 'bold');

if MakeVideo
    set(gcf, 'Position', [300, 300, 700, 600]);
    drawnow;
    
    if exist(VideoName, 'file')
        delete(VideoName);
    end
    
    v = VideoWriter(VideoName, 'MPEG-4');
    v.FrameRate = fps;
    open(v);
    
    idx0 = frame_indices(1);
    set(h_rod,  'XData', x0_vals(idx0,:), 'YData', z0_vals(idx0,:));
    set(h_root, 'XData', x0_vals(idx0,1), 'YData', z0_vals(idx0,1));
    set(h_tip,  'XData', x0_vals(idx0,end), 'YData', z0_vals(idx0,end));
    set(h_time, 'String', sprintf('Time: %.2f s', t(idx0)));
    drawnow;
    writeVideo(v, getframe(gcf));
    
    fprintf('正在生成视频文件 %s ...\n', VideoName);

    for k = 2:length(frame_indices)
        idx = frame_indices(k);
        set(h_rod,  'XData', x0_vals(idx,:), 'YData', z0_vals(idx,:));
        set(h_root, 'XData', x0_vals(idx,1), 'YData', z0_vals(idx,1));
        set(h_tip,  'XData', x0_vals(idx,end), 'YData', z0_vals(idx,end));
        set(h_time, 'String', sprintf('Time: %.2f s', t(idx)));
        drawnow;
        writeVideo(v, getframe(gcf));
    end
    
    close(v);
    fprintf('视频生成完毕！保存至: %s\n', VideoName);
else
    idx0 = frame_indices(1);
    set(h_rod,  'XData', x0_vals(idx0,:), 'YData', z0_vals(idx0,:));
    set(h_root, 'XData', x0_vals(idx0,1), 'YData', z0_vals(idx0,1));
    set(h_tip,  'XData', x0_vals(idx0,end), 'YData', z0_vals(idx0,end));
    set(h_time, 'String', sprintf('Time: %.2f s', t(idx0)));
    drawnow;
end

fprintf('========== 仿真完美结束 ==========\n');
end

%% ========== 子函数：支持全矩阵计算的 O(n) 积分 ==========
function I = fast_double_integral(w, C, v)
    wv = w .* v;
    wCv = (w .* C) .* v;
    prefix = cumsum(wv, 1);
    suffix = cumsum(wCv, 1, 'reverse');
    left = C .* (prefix - wv);
    I = left + suffix;
end

%% ========== 子函数 1: PDE 求解器 (向量化 + 弱态质量矩阵) ==========
function [t, Y_full, S] = solvePDE_MassMatrix_Vectorized(h0, eps, nGrid, Tmax, F, beta)
    L = 1; ds = L / nGrid; S = linspace(0, L, nGrid+1)'; nNodes = length(S);
    
    w = ones(nNodes,1);
    if mod(nGrid,2) == 0; w(2:2:end-1) = 4; w(3:2:end-2) = 2; w = w * ds / 3;
    else; w(1) = ds/2; w(end) = ds/2; end
    C = cumsum(w, 'reverse');
    
    Wdouble = zeros(nNodes, nNodes);
    for i = 1:nNodes; for j = 1:nNodes; Wdouble(i,j) = w(j) * sum(w(max(i,j):end)); end; end
    
    D1 = sparse(nNodes, nNodes); D1(1, 1:5) = [-25, 48, -36, 16, -3] / (12*ds); D1(2, 1:5) = [-3, -10, 18, -6, 1] / (12*ds);
    for i = 3:nNodes-2; D1(i, i-2:i+2) = [1, -8, 0, 8, -1] / (12*ds); end; 
    D1(nNodes-1, nNodes-4:nNodes) = [-1, 6, -18, 10, 3] / (12*ds);
    
    D2 = sparse(nNodes, nNodes); D2(2, 1:6) = [10, -15, -4, 14, -6, 1] / (12*ds^2);
    for i = 3:nNodes-2; D2(i, i-2:i+2) = [-1, 16, -30, 16, -1] / (12*ds^2); end
    D2(nNodes-1, nNodes-4:nNodes) = [-1, 4, 6, -20, 11] / (12*ds^2); 
    D2(nNodes, nNodes-2:nNodes) = [-2, 32, -30] / (12*ds^2);
    
    free = 2:nNodes; nFree = length(free); y0_free = zeros(2*nFree, 1);
    I_reg = 1e-9 * speye(nFree); 

    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-5, ...
                     'Mass', @mass_matrix, 'MStateDependence', 'weak', ...
                     'Vectorized', 'on', 'MaxOrder', 2); 
    
    [t, Y_free] = ode15s(@odeRHS, [0, Tmax], y0_free, options);
    
    nt = length(t); Y_full = zeros(nt, 2*nNodes);
    Y_full(:, free) = Y_free(:, 1:nFree);             
    Y_full(:, nNodes + free) = Y_free(:, nFree+1:end); 

    function res = odeRHS(t, y_free)
        K = size(y_free, 2);
        g = zeros(nNodes, K); g(free, :) = y_free(1:nFree, :);
        v = zeros(nNodes, K); v(free, :) = y_free(nFree+1:end, :);
        cosg = cos(g); sing = sin(g);
        
        Nl10    = F.Nl10(S, t);    if isscalar(Nl10), Nl10 = Nl10 * ones(nNodes, 1); end
        Nl10_S  = F.Nl10_S(S, t);  if isscalar(Nl10_S), Nl10_S = Nl10_S * ones(nNodes, 1); end
        Nl10_t  = F.Nl10_t(S, t);  if isscalar(Nl10_t), Nl10_t = Nl10_t * ones(nNodes, 1); end
        Nl10_tt = F.Nl10_tt(S, t); if isscalar(Nl10_tt), Nl10_tt = Nl10_tt * ones(nNodes, 1); end
        Nl11    = F.Nl11(S, t);    if isscalar(Nl11), Nl11 = Nl11 * ones(nNodes, 1); end
        Nl11_S  = F.Nl11_S(S, t);  if isscalar(Nl11_S), Nl11_S = Nl11_S * ones(nNodes, 1); end
        Nl11_t  = F.Nl11_t(S, t);  if isscalar(Nl11_t), Nl11_t = Nl11_t * ones(nNodes, 1); end
        
        gSS = D2 * g; gS = D1 * g; 
        vSS = D2 * v; vS = D1 * v; 

        B_val   = -Nl11(end);
        B_val_t = -Nl11_t(end);

        gS(end, :) = B_val;
        vS(end, :) = B_val_t;

        gSS(end, :) = (2 * g(end-1, :) - 2 * g(end, :) + 2 * ds * B_val) / ds^2;
        vSS(end, :) = (2 * v(end-1, :) - 2 * v(end, :) + 2 * ds * B_val_t) / ds^2;
        
        gSS_eff = gSS + beta * vSS;
        gS_eff  = gS  + beta * vS;
        
        PartA = ( 8*h0^2 .* ( -(Nl11 + gS_eff).*Nl10_S + Nl10.*(Nl11_S + gSS_eff) ) ) ./ ( 3 * Nl10.^3 );
        
        term_Ix_v = cosg .* (Nl10_tt - Nl10 .* v.^2) - sing .* (2.*v.*Nl10_t);
        term_Iz_v = sing .* (Nl10_tt - Nl10 .* v.^2) + cosg .* (2.*v.*Nl10_t);
        
        Ix_v = fast_double_integral(w, C, term_Ix_v); 
        Iz_v = fast_double_integral(w, C, term_Iz_v);
        
        R = -PartA - eps*sing.*(Ix_v + (1-S)*F.cxddot(t)) + eps*cosg.*(Iz_v + (1-S)*F.czddot(t));
        res = [v(free, :); -R(free, :)];
    end

    function M_sys = mass_matrix(t, y_free)
        g = zeros(nNodes, 1); g(free) = y_free(1:nFree);
        Nl10 = F.Nl10(S, t); if isscalar(Nl10), Nl10 = Nl10 * ones(nNodes, 1); end
        
        M_full_Positive = eps * Wdouble .* (Nl10' .* cos(g - g'));
        M_sub = M_full_Positive(free, free) + I_reg;
        M_sys = blkdiag(speye(nFree), M_sub);
    end
end

%% ========== 子函数 2: 后处理 (顺序代入) ==========
function [g_tt, g_S, f_vals, x0_vals, z0_vals] = compute_kinematics_PostProcess(t, S, g_sol, g_t, h0, eps, F, PlotMode)
    nt = length(t); ns = length(S); ds = S(2) - S(1);
    w = ones(ns,1);
    if mod(ns-1, 2) == 0; w(2:2:end-1) = 4; w(3:2:end-2) = 2; w = w * ds / 3;
    else; w(1) = ds/2; w(end) = ds/2; end
    C = cumsum(w, 'reverse');
    Wdouble = zeros(ns, ns);
    for i = 1:ns; for j = 1:ns; Wdouble(i,j) = w(j) * sum(w(max(i,j):end)); end; end
    
    D1 = sparse(ns, ns); D1(1, 1:5) = [-25, 48, -36, 16, -3] / (12*ds); D1(2, 1:5) = [-3, -10, 18, -6, 1] / (12*ds);
    for i = 3:ns-2; D1(i, i-2:i+2) = [1, -8, 0, 8, -1] / (12*ds); end; 
    D1(ns-1, ns-4:ns) = [-1, 6, -18, 10, 3] / (12*ds);
    
    D2 = sparse(ns, ns); D2(2, 1:6) = [10, -15, -4, 14, -6, 1] / (12*ds^2);
    for i = 3:ns-2; D2(i, i-2:i+2) = [-1, 16, -30, 16, -1] / (12*ds^2); end
    D2(ns-1, ns-4:ns) = [-1, 4, 6, -20, 11] / (12*ds^2); 
    D2(ns, ns-2:ns) = [-2, 32, -30] / (12*ds^2);
    
    g_tt = zeros(nt, ns); g_S  = zeros(nt, ns); f_vals = zeros(nt, ns); x0_vals = zeros(nt, ns); z0_vals = zeros(nt, ns);
    free = 2:ns; I_reg = 1e-9 * speye(length(free));
    
    for i = 1:nt
        ti = t(i); g = g_sol(i,:)'; v = g_t(i,:)'; cosg = cos(g); sing = sin(g);
        
        Nl10 = F.Nl10(S, ti); if isscalar(Nl10), Nl10 = Nl10 * ones(ns, 1); end
        Nl10_S = F.Nl10_S(S, ti); if isscalar(Nl10_S), Nl10_S = Nl10_S * ones(ns, 1); end
        Nl10_t = F.Nl10_t(S, ti); if isscalar(Nl10_t), Nl10_t = Nl10_t * ones(ns, 1); end
        Nl10_tt = F.Nl10_tt(S, ti); if isscalar(Nl10_tt), Nl10_tt = Nl10_tt * ones(ns, 1); end
        Nl11 = F.Nl11(S, ti); if isscalar(Nl11), Nl11 = Nl11 * ones(ns, 1); end
        Nl11_S = F.Nl11_S(S, ti); if isscalar(Nl11_S), Nl11_S = Nl11_S * ones(ns, 1); end
        Nl12 = F.Nl12(S, ti); if isscalar(Nl12), Nl12 = Nl12 * ones(ns, 1); end
        
        gSS = D2 * g; gS = D1 * g; 

        B_val = -Nl11(end);
        gS(end) = B_val; 
        gSS(end) = (2 * g(end-1) - 2 * g(end) + 2 * ds * B_val) / ds^2;
        g_S(i,:) = gS';
        
        PartA = ( 8*h0^2 .* ( -(Nl11 + gS).*Nl10_S + Nl10.*(Nl11_S + gSS) ) ) ./ ( 3 * Nl10.^3 );
        term_Ix_v = cosg .* (Nl10_tt - Nl10 .* v.^2) - sing .* (2.*v.*Nl10_t);
        term_Iz_v = sing .* (Nl10_tt - Nl10 .* v.^2) + cosg .* (2.*v.*Nl10_t);
        
        Ix_v = fast_double_integral(w, C, term_Ix_v); Iz_v = fast_double_integral(w, C, term_Iz_v);
        R = -PartA - eps*sing.*(Ix_v + (1-S)*F.cxddot(ti)) + eps*cosg.*(Iz_v + (1-S)*F.czddot(ti));
        
        M_full_Positive = eps * Wdouble .* (Nl10' .* cos(g - g'));
        a_free = (M_full_Positive(free, free) + I_reg) \ (-R(free));
        a = zeros(ns,1); a(free) = a_free; g_tt(i,:) = a';
        
        Ix_full = cosg .* (Nl10_tt - Nl10 .* v.^2) - sing .* (2.*v.*Nl10_t + Nl10.*a);
        Iz_full = sing .* (Nl10_tt - Nl10 .* v.^2) + cosg .* (2.*v.*Nl10_t + Nl10.*a);      
        Ix_double = fast_double_integral(w, C, Ix_full); Iz_double = fast_double_integral(w, C, Iz_full);
        
        Part1 = 3 * cosg .* ( - (Ix_double + (1-S).*F.cxddot(ti)) );
        Part2 = 3 * sing .* ( - (Iz_double + (1-S).*F.czddot(ti)) );
        B_term = 2*h0 * (Nl11 + gS) .* (Nl11 + 3*gS) + Nl10 .* (3*Nl11 + 2*h0*Nl12 + 3*gS);
        Part3  = (8 * h0 * B_term) ./ (eps * Nl10.^2);
        
        f = Nl10 .* (1 + eps * (1/24) * (Part1 + Part2 + Part3)); f_vals(i,:) = f';
        fx = f .* cosg; fz = f .* sing; rel_x = cumtrapz(S, fx)'; rel_z = cumtrapz(S, fz)';
        
        if strcmp(PlotMode, 'Absolute'); x0_vals(i,:) = rel_x + F.cx(ti); z0_vals(i,:) = rel_z + F.cz(ti);
        else; x0_vals(i,:) = rel_x; z0_vals(i,:) = rel_z; end
    end
end