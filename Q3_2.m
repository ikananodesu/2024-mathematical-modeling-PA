clc, clear, close all;
% 输入运行时间t
flag=1;
P=45.03;%最小间距
    figure(1);
    hold on;
theta=20*2*pi;
% 绘制阿基米德螺线
theta_vals = linspace(0, max(theta), 1000); % 设置theta范围
R_vals = get_length(theta_vals,P); % 根据theta计算螺线的半径
X_vals = R_vals .* cos(theta_vals); % 计算螺线的X坐标
Y_vals = R_vals .* sin(theta_vals); % 计算螺线的Y坐标
plot(X_vals, Y_vals, 'k-', 'LineWidth', 1); % 绘制螺线，使用黑色虚线
t=221.5;%结束时间
while(flag)
    [L, A, X, Y, V] = Q1_f(P,t);
    X = X';
    Y = Y';
    V = V';
    % 绘制螺线点图
    scatter(X, Y, 'filled');
        % 画半径为450的圆
    theta_circle = linspace(0, 2 * pi, 5000); % 0到2pi的角度
    r_circle = 450; % 圆的半径
    x_circle = r_circle * cos(theta_circle); % 圆上x坐标
    y_circle = r_circle * sin(theta_circle); % 圆上y坐标
    plot(x_circle, y_circle, 'm-', 'LineWidth', 2); % 用蓝色虚线绘制圆
    axis equal;
    %初始化矩形存储
    rectangles = [];
    % 绘制相邻点之间的矩形
    for i = 1:length(X) - 1
        % 计算相邻把手之间的向量
        dx = X(i + 1) - X(i);
        dy = Y(i + 1) - Y(i);
        length_vec = sqrt(dx^2 + dy^2);

        % 计算向量的单位方向
        unit_vec = [dx, dy] / length_vec;

        % 垂直向量的单位方向（旋转90度）
        perp_vec = [-unit_vec(2), unit_vec(1)];

        % 定义矩形的四个顶点
        p1 = [X(i), Y(i)] + perp_vec * 15 - unit_vec * 27.5;  % 上边
        p2 = [X(i + 1), Y(i + 1)] + perp_vec * 15 + unit_vec * 27.5;  % 上边
        p3 = [X(i + 1), Y(i + 1)] - perp_vec * 15 + unit_vec * 27.5;  % 下边
        p4 = [X(i), Y(i)] - perp_vec * 15 - unit_vec * 27.5;  % 下边

        % 将顶点坐标存入矩形数组
        rectangles = [rectangles; struct('points', [p1; p2; p3; p4], 'index', i)];

        % 绘制矩形
         fill([p1(1), p2(1), p3(1), p4(1)], [p1(2), p2(2), p3(2), p4(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
    end
    % 碰撞检测：检查非相邻矩形是否重叠
    for i = 1:length(rectangles)
        for j = i + 2:length(rectangles)  % 检查非相邻的矩形
            if checkOverlap(rectangles(i).points, rectangles(j).points)
                disp(['检测到矩形 ', num2str(rectangles(i).index), ' 和 ', num2str(rectangles(j).index), ' 发生重叠！']);
                % 标注重叠矩形
                fill(rectangles(i).points(:, 1), rectangles(i).points(:, 2), 'g', 'FaceAlpha', 0.4);
                fill(rectangles(j).points(:, 1), rectangles(j).points(:, 2), 'g', 'FaceAlpha', 0.4);
                flag=0;
            end
        end
    end
    scatter(0, 0, 'filled','r');
    scatter(X(2:end), Y(2:end),10, 'filled', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    scatter(X(1), Y(1),10, 'filled', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    grid on
    hold off;
    if(X(1)^2+Y(1)^2<=450^2)
        flag=0;
    end
    t=t+0.5;
end
    xlabel('X 坐标/cm');
    ylabel('Y 坐标/cm');
%% 极限圆进行碰撞检测
function [L, A, X, Y, V] = Q1_f(P,t)
    % 定义阿基米德螺线和龙头/龙身的计算方法
    a = 0; % 起始半径
    b = P / (2 * pi); % 间距系数
    body_num = 223; % 板凳总量
    theta = zeros(1, body_num + 1);
    V = 100; % 龙头的速度是100cm/s
    S = V * t; % 龙头运行距离
    theta_start = 16 * 2 * pi; % 出发点是第16圈
    theta(1) = solve_theta_from_arc_length(S, theta_start,P); % 初始化龙头点的位置
    lp_head = 286; % 龙头的板长
    lp_body = 165; % 龙身的板长
    % 计算每个点的位置
    for i = 1:body_num
        if i == 1
            theta(i + 1) = solve_theta2(theta(i), lp_head,P);
        else
            theta(i + 1) = solve_theta2(theta(i), lp_body,P);
        end
    end
    R = get_length(theta,P);
    V = zeros(1, body_num + 1);
    V(1) = 100;
    for i = 1:body_num
        V(i + 1) = V(i) * R(i) /R(i + 1);
    end
    X = R .* cos(theta);
    Y = R .* sin(theta);
    for i = 2:5
        L(i - 1) = sqrt((X(1) - X(i))^2 + (Y(1) - Y(i))^2);
    end
    lp_head = 286;
    lp_body = 165;
    beta1 = atand(15 / lp_head);
    beta2 = atand(15 / lp_body);
    AB = [X(1) - X(2), Y(1) - Y(2)];
    BC = [X(2) - X(3), Y(2) - Y(3)];
    A = dot(AB, BC);
    mAB = norm(AB);
    mBC = norm(BC);
    Angle = A / (mAB * mBC);
    Angle = acosd(Angle);
    if (beta1 + beta2 >= Angle)
        A = true;
    else
        A = false;
    end
end

% 碰撞检测函数，检查两个矩形是否重叠
function overlap = checkOverlap(points1, points2)
    % 使用分离轴定理（SAT）来判断两个旋转矩形是否重叠
    overlap = false;
    % 获取所有边作为分离轴
    edges = [points1(2,:) - points1(1,:); points1(3,:) - points1(2,:); points2(2,:) - points2(1,:); points2(3,:) - points2(2,:)];
    for i = 1:size(edges, 1)
        % 计算分离轴
        axis = [-edges(i,2), edges(i,1)];
        % 归一化分离轴
        axis = axis / norm(axis);
        % 在分离轴上投影所有点
        proj1 = points1 * axis';
        proj2 = points2 * axis';
        % 检查投影是否重叠
        if max(proj1) < min(proj2) || max(proj2) < min(proj1)
            return; % 如果发现无重叠，则没有碰撞
        end
    end
    overlap = true; % 如果所有投影都有重叠，说明有碰撞
end

% 其他辅助函数保留
function theta2 = solve_theta2(theta1, X,P)
    b = P / (2 * pi); % 间距系数
    L1 = get_length(theta1,P);
    theta2_initial_guess = theta1 + X / L1; % 初始猜测值
    L2 = @(theta2) get_length(theta2,P);
    distance_function = @(theta2) sqrt(L1^2 + L2(theta2)^2 ...
        - 2 * L1 * L2(theta2) * cos(theta1 - theta2)) - X;
    options = optimoptions('fsolve', 'Display', 'off');
    theta2 = fsolve(distance_function, theta2_initial_guess, options);
    if theta2 <= theta1
        warning('求解得到的theta2 <= theta1，调整初始猜测值后重试');
        theta2_initial_guess = theta1 + abs(theta2_initial_guess - theta1) + 1e-6; % 增大初始猜测值
        theta2 = fsolve(distance_function, theta2_initial_guess, options);
    end
end

function R = get_length(theta,P)
    a = 0; % 起始半径
    b = P / (2 * pi); % 间距系数
    R = a + b * theta; % 螺线的半径公式
end

function theta = solve_theta_from_arc_length(L, theta_start,P)
    a = 0; % 起始半径
    b = P / (2 * pi); % 间距系数
    arc_length_diff = @(theta) get_arc_length(theta_start, theta,P) - L;
    theta_initial_guess = theta_start - L / b;
    options = optimoptions('fsolve', 'Display', 'off');
    theta = fsolve(arc_length_diff, theta_initial_guess, options);
end

function L = get_arc_length(theta_start, theta,P)
    a = 0; % 起始半径
    b = P / (2 * pi); % 间距系数
    arc_length = @(theta) sqrt(b^2 + (a + b * theta).^2);
    L = integral(arc_length, theta, theta_start);
end