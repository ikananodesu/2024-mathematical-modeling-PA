clc; clear; close all;

% 设定圆的半径
radius = 450;
theta_circle = linspace(0.81*pi, 10*pi, 10000); % 生成从 0 到 2*pi 的点

% 计算圆的坐标
x_circle = radius * cos(theta_circle);
y_circle = radius * sin(theta_circle);

% 绘制圆
figure;
plot(x_circle, y_circle, 'b-', 'LineWidth', 2);
hold on;

% 设定阿基米德螺线的参数
L = 3000;
a = 0;
b = 170;
theta_spiral = linspace(0.81*pi, 10*pi, 10000); % 螺旋线旋转了4圈
r_spiral = a + b * theta_spiral;
x_spiral = r_spiral .* cos(theta_spiral);
y_spiral = r_spiral .* sin(theta_spiral);
% 绘制阿基米德螺线
plot(x_spiral, y_spiral, 'r-', 'LineWidth', 2);

% 绘制反向阿基米德螺线
x_spiral_reverse = -x_spiral;
y_spiral_reverse = -y_spiral;
plot(x_spiral_reverse, y_spiral_reverse, 'g-', 'LineWidth', 2);

% 求交点
theta_intersect_pos = radius / b;  % 交点处的theta
r_intersect_pos = a + b * theta_intersect_pos;
x_intersect_pos = r_intersect_pos * cos(theta_intersect_pos);
y_intersect_pos = r_intersect_pos * sin(theta_intersect_pos);

% 计算法向量和延长线
dx_dtheta_pos = -r_intersect_pos * sin(theta_intersect_pos) + b * cos(theta_intersect_pos);
dy_dtheta_pos = r_intersect_pos * cos(theta_intersect_pos) + b * sin(theta_intersect_pos);
normal_vector_pos = [-dy_dtheta_pos, dx_dtheta_pos]; % 垂直于切向量
normal_vector_pos = normal_vector_pos / norm(normal_vector_pos);  % 单位化法向量

% 求解法向量延长线与圆的交点
syms t;
x_line = x_intersect_pos + t * normal_vector_pos(1);
y_line = y_intersect_pos + t * normal_vector_pos(2);
eqn = (x_line)^2 + (y_line)^2 == radius^2;
sol_t = double(solve(eqn, t));
t_valid_pos = max(sol_t);

% 绘制法向量延长线上的点
x_end_pos = x_intersect_pos + t_valid_pos * normal_vector_pos(1);
y_end_pos = y_intersect_pos + t_valid_pos * normal_vector_pos(2);
t_values = linspace(0, t_valid_pos, L);
x_points = x_intersect_pos + t_values * normal_vector_pos(1);
y_points = y_intersect_pos + t_values * normal_vector_pos(2);

% 反向螺旋线法向量
x_intersect_neg = -x_intersect_pos;
y_intersect_neg = -y_intersect_pos;
normal_vector_neg = -normal_vector_pos;

% 求解反向法向量延长线与圆的交点
x_line_neg = x_intersect_neg + t * normal_vector_neg(1);
y_line_neg = y_intersect_neg + t * normal_vector_neg(2);
eqn_neg = (x_line_neg)^2 + (y_line_neg)^2 == radius^2;
sol_t_neg = double(solve(eqn_neg, t));
t_valid_neg = max(sol_t_neg);

% 绘制反向法向量延长线上的点
x_end_neg = x_intersect_neg + t_valid_neg * normal_vector_neg(1);
y_end_neg = y_intersect_neg + t_valid_neg * normal_vector_neg(2);
t_values_neg = linspace(0, t_valid_neg, L);
x_points_neg = x_intersect_neg + t_values_neg * normal_vector_neg(1);
y_points_neg = y_intersect_neg + t_values_neg * normal_vector_neg(2);

%% 计算 R1 和 R2 及其圆心位置
left_point = [x_points; y_points];
right_point = [x_points_neg; y_points_neg];
l_begin_in = left_point(:, 1);
r_begin_in = right_point(:, 1);
l_begin_end = left_point(:, end);
r_begin_end = right_point(:, end);

% R1 半径的计算
R1 = zeros(1, L);
for i = 1:L
    R1(i) = sqrt((left_point(1, i) - l_begin_in(1))^2 + (left_point(2, i) - l_begin_in(2))^2);
end

% 计算左点和右点之间的距离矩阵
S = zeros(L, L);
R2 = zeros(1, L);
for i = 1:L
    for j = 1:L
        S(i,j) = sqrt((left_point(1, i) - right_point(1, j))^2 + (left_point(2, i) - right_point(2, j))^2);
        if abs(R1(i) + R1(j) - S(i,j)) <= 1
            label(i, j) = 1;
        end
    end
end

% 找到最接近的点对
[i, j] = find(label == 1);  % 找到满足条件的点
left_point_selected = left_point(:, i);
right_point_selected = right_point(:, j);
R2 = R1(j);
R1 = R1(i);

% 计算两圆心的平均位置，并求得最小的点
sum_ave = left_point_selected + right_point_selected;
S = sum_ave.^2;
SS = sum(S,1);
[A] = find(SS == min(SS));
left_point_place = left_point_selected(:,mean(A));
right_point_place = right_point_selected(:,mean(A));

% 绘制求解得到的两个圆心和其圆
R_avg = R1(mean(A));  % 使用 R1 作为两个圆的半径

% 计算 x_intersect_pos 和 y_intersect_pos 在左圆上对应的角度
theta_pos = atan2(y_intersect_pos - left_point_place(2), x_intersect_pos - left_point_place(1));
theta_pos_deg = rad2deg(theta_pos); % 转换为角度
theta_posA = atan2(left_point_place(2),left_point_place(1));
theta_posA_deg = rad2deg(theta_posA); % 转换为角度
theta = linspace(theta_pos_deg*pi/180, 0.9*pi-theta_posA_deg*pi/180, 2000);

% 输出角度
disp(['交点在左圆上的角度为：', num2str(theta_pos_deg), ' 度']);

% 绘制以 left_point_place 为圆心的圆
x_left_circle = R_avg * cos(theta) + left_point_place(1);
y_left_circle = R_avg * sin(theta) + left_point_place(2);
plot(x_left_circle, y_left_circle, 'm-', 'LineWidth', 2);

X=(x_intersect_neg+x_intersect_pos)/2;
Y=(y_intersect_neg+y_intersect_pos)/2;

RX=2*X-x_left_circle;
RY=2*Y-y_left_circle;
RX=flip(RX);
RY=flip(RY);
plot(RX, RY, 'm-', 'LineWidth', 2);
% 设置图形参数
axis equal;
title('螺旋线、圆及计算所得圆心');
xlabel('X轴');
ylabel('Y轴');
legend('圆', '正阿基米德螺线', '反阿基米德螺线', '交点', '法向量散点', 'Left Circle', 'Right Circle');
hold off;

XR=[x_left_circle(1:1987),(x_left_circle(1989:end)+RX(1:12))/2,RX(14:end)];
YR=[y_left_circle(1:1987),(y_left_circle(1989:end)+RY(1:12))/2,RY(14:end)];

x_spiral_reverse = -x_spiral;
y_spiral_reverse = -y_spiral;
%10000 3986 10000
% ans_get_x=[flip(x_spiral),XR,x_spiral_reverse];
% ans_get_y=[flip(y_spiral),YR,y_spiral_reverse];
ans_get_x=[flip(x_spiral),XR,x_spiral_reverse];
ans_get_y=[flip(y_spiral),YR,y_spiral_reverse];
way=[ans_get_x;ans_get_y];
plot(ans_get_x,ans_get_y)
save('问题四路径','way');