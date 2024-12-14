clc, clear, close all;
% 加载路径数据
way = load("问题四路径.mat");
way = way.way; % 'way' 包含路径的 X 和 Y 坐标

% 初始化参数
target_X = -396.085526860245; % 初始X坐标
target_Y = 213.579623119440; % 初始Y坐标
r = 241.1515; % 圆半径
head = 286;
body = 165;
body_num = 223; % 板凳总量
points = zeros(2, body_num+1); % 储存身体各部分坐标
label = zeros(1, body_num+1); % 路径映射标签
len_points = zeros(1, body_num+1); % 距离出发点的路径长度
num_points = size(way, 2); % 地图标签路径区间
V = zeros(1, body_num+1); % 速度
a = zeros(1, body_num+1); % 角速度
V(1) = 100; % 龙头初始速度

% 设定时间
t_min = -100; % 最小时间
t_max = 100; % 最大时间
dt = 1; % 时间步长
timesteps = t_min:dt:t_max; % 时间步
S = 0; % 初始线路程

% 用于保存每个时刻的points和速度V
points_over_time = zeros(2 * (body_num + 1), length(timesteps)); 
V_over_time = zeros(body_num + 1, length(timesteps)); 

% 计算出发点到路径的最近点
dist_to_target = sqrt((way(1, :) - target_X).^2 + (way(2, :) - target_Y).^2);
[~, begin_index] = min(dist_to_target); % 最近的路径标签

for t_idx = 1:length(timesteps)
    t = timesteps(t_idx);    
    % 计算当前时间的行程
    S = V(1) * t; % 运行距离
    [lengthA, i] = func_get_len(S, begin_index, way); % 出发点-龙头位置的弧积分
    points(:, 1) = way(:, i); % 存储龙头坐标
    len_points(1) = lengthA;
    label(1) = i;
    
    %% 回溯求解身体坐标
    for j = 1:body_num
        if j == 1
            L = head; % 龙头部分
        else
            L = body; % 其余部分
        end
        dir = -1; % 回溯方向
        [points(1, j+1), points(2, j+1), label(j+1)] = select_f(points(1, j), points(2, j), label(j), L, dir, way);
    end
    
    % 将 points 的 [X, Y] 坐标平展为 [X1; Y1; X2; Y2; ...]
    flat_points = reshape(points, [], 1);
    
    % 保存每个时间步的平展后的 points
    points_over_time(:, t_idx) = flat_points;
    
    %% 分区规划计算速度
    for k = 1:body_num
        if label(k) <= 10000 % 入旋线
            V(k+1) = sqrt(points(1,k+1)^2 + points(2,k+1)^2) / sqrt(points(1,k)^2 + points(2,k)^2) * V(k);
            a(k+1) = V(k+1) / sqrt(points(1,k+1)^2 + points(2,k+1)^2);
        elseif label(k) <= 13986 && label(k) > 10000 % 圆弧部分
            V(k+1) = V(k);
            a(k+1) = V(k+1) / r;
        else % 出旋线
            V(k+1) = sqrt(points(1,k+1)^2 + points(2,k+1)^2) / sqrt(points(1,k)^2 + points(2,k)^2) * V(k);
            a(k+1) = V(k+1) / sqrt(points(1,k+1)^2 + points(2,k+1)^2);
        end
    end
    
    % 保存速度
    V_over_time(:, t_idx) = V;
end

%% 绘制路径和长方形板凳
figure;
hold on;
theta_circle = linspace(0, 2*pi, 2000);
r_circle = 450; % 大圆半径
x_circleA = r_circle * cos(theta_circle);
y_circleA = r_circle * sin(theta_circle);
plot(x_circleA, y_circleA, 'm-', 'LineWidth', 1.5);
plot(way(1,1:10000), way(2,1:10000), 'r', 'LineWidth', 1.5);
plot(way(1,10001:13988), way(2,10001:13988), 'c-', 'LineWidth', 1.5);
plot(way(1,13989:end), way(2,13989:end), 'g', 'LineWidth', 1.5);

% 绘制身体坐标点
scatter(points(1,2:end), points(2,2:end), 20, 'filled', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
plot(points(1,:), points(2,:), 'b.-');
axis equal;
grid on;
xlabel('X坐标/cm');
ylabel('Y坐标/cm');
hold off;

figure;
surf(V_over_time);
shading interp
%% 函数定义部分
% 选择并回溯查找坐标
function [X,Y,label_out] = select_f(points_x,points_y,label_in,L,dir,way)
    i = label_in;
    if dir == 1 % 向前搜索
        while i < length(way) && abs(sqrt((way(1,i) - points_x)^2 + (way(2,i) - points_y)^2) - L) > 10
            label_out = i;
            i = i + 1;
        end
    else % 向后搜索
        while i > 0 && abs(sqrt((way(1,i) - points_x)^2 + (way(2,i) - points_y)^2) - L) > 10
            label_out = i;
            i = i - 1;
        end
    end
    X = way(1,label_out);
    Y = way(2,label_out);
end

% 曲线积分反解路径点
function [length, i] = func_get_len(S, begin_index, way)
    num_points = size(way, 2);
    length = 0;
    i = begin_index; % 出发点坐标
    if S > 0 % 积分大于0，往前找
        while i < num_points && length < S
            dX = way(1,i+1) - way(1,i);
            dY = way(2,i+1) - way(2,i);
            length = length + sqrt(dX^2 + dY^2);
            i = i + 1;
        end
    elseif S < 0 % 积分小于0，往后找
        while i > 1 && length < abs(S)
            dX = way(1,i-1) - way(1,i);
            dY = way(2,i-1) - way(2,i);
            length = length + sqrt(dX^2 + dY^2);
            i = i - 1;
        end
    end
end


