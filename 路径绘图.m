clc, clear, close all;
way = load("问题四路径.mat");
way = way.way;
%% 绘制不同区域的路径
figure(1);
hold on;
plot(way(1, 1:10000), way(2, 1:10000), 'r-', 'LineWidth', 2.5); % 盘入阶段
plot(way(1, 10001:11993), way(2, 10001:11993), 'c-', 'LineWidth', 2.5); % 圆弧1阶段
plot(way(1, 11994:13986), way(2, 11994:13986), 'b-', 'LineWidth', 2.5); % 圆弧2阶段
plot(way(1, 13987:end), way(2, 13987:end), 'g-', 'LineWidth', 2.5); % 盘出阶段
plot(way(1, 10000), way(2, 10000), 'ko', 'MarkerSize', 3, 'LineWidth', 2, 'MarkerFaceColor', 'k');
grid on;
legend('盘入阶段', '圆弧1阶段','圆弧2阶段','盘出阶段','零时刻点');
xlabel('X 坐标/cm');
ylabel('Y 坐标/cm');
axis equal;
hold off;

surf(V)
shading interp
xlabel('运行时间t');
ylabel('板凳把手坐标i');
zlabel('速度m/s');

data=xlsread("C:\Users\86188\Desktop\24国赛数模 伊卡\NO33\Q4数据.xlsx")
data=data/100;
figure;
hold on;
grid on;
plot(data(:,1),'b-','LineWidth',2);
plot(data(:,2),'r-','LineWidth',2);
plot(data(:,3),'g-','LineWidth',2);
xlabel('关节编号i');
ylabel('速度v m/s');
legend('t=-10 盘入阶段','t=10 调头阶段','t=100 盘出阶段')
figure;
hold on;
grid on;
plot(diff(data(:,1)),'b-','LineWidth',2);
plot(diff(data(:,2)),'r-','LineWidth',2);
plot(diff(data(:,3)),'g-','LineWidth',2);
xlabel('关节编号i');
ylabel('加速度a m/s2');
legend('t=-10 盘入阶段','t=10 调头阶段','t=100 盘出阶段')
