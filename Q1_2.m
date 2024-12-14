clc, clear, close all;
place2=[];
Vans=[];
%搜索时刻-把手-速度
for t=0:300
    [V,X,Y]=func(t);
    place1=[];
    for i=1:length(X)
        place1=[place1;X(i);Y(i)];
    end
    place2=[place2,place1/100];
    Vans=[Vans,V/100];
end
figure;
surf(Vans);
xlabel('运行时间t');
ylabel('把手编号');
zlabel('速度m/s');
grid on;
shading interp;

function [V,X,Y] = func(t)
P=55;%螺距55cm
a=0;
b=P/(2*pi);%间距系数
body_num=223;%板凳总量
theta=zeros(1,body_num+1);%初始化极坐标存储
V=100;%龙头速度100cm/s
S=V*t;%龙头移动弧长
theta_start=16*2*pi;%第16圈出发
theta(1)=solve_theta_from_arc_length(S,theta_start);%初始化龙头位置
lp_head=286;%板子参数
lp_body=165;
%迭代把手角坐标
for i=1:body_num
    if i==1
        theta(i+1)=solve_theta2(theta(i),lp_head);
    else
        theta(i+1)=solve_theta2(theta(i),lp_body);
    end
end
%通过阿基米德螺旋线反解与原点的半径
R=get_length(theta);
V=zeros(1, body_num + 1);%初始化存储
V(1)=100;%龙头速度
%迭代把手速度
for i=1:body_num
    V(i+1)=V(i)*cos((theta(i+1)-theta(i))/180);
end
V=V';
X=R.*cos(theta);
Y=R.*sin(theta);
X=X';
Y=Y';
end
% 计算尾部节点角度的函数
function theta2 = solve_theta2(theta1, X)
P=55;%螺距55cm
a=0;
b=P/(2*pi);%间距系数
L1=get_length(theta1);
%搜索拟合角度
theta2_initial_guess=theta1+X/L1;
L2 = @(theta2) get_length(theta2);
distance_function=@(theta2) sqrt(L1^2 +L2(theta2)^2-2*L1*L2(theta2)*cos(theta1-theta2))-X;
options=optimoptions('fsolve','Display','off');
theta2=fsolve(distance_function,theta2_initial_guess,options);
if theta2<=theta1
    theta2_initial_guess=theta1+abs(theta2_initial_guess-theta1)+1e-6;%增加步长
    theta2=fsolve(distance_function, theta2_initial_guess, options);
end
end

%计算螺旋线半径
function R=get_length(theta)
P=55;%螺距55cm
a=0;
b=P/(2*pi);%间距系数
R=a+b*theta;
end

% 弧长反解转角
function theta = solve_theta_from_arc_length(L,theta_start)
P=55;%螺距55cm
a=0;
b=P/(2*pi);%间距系数
arc_length_diff=@(theta) get_arc_length(theta_start,theta)-L;
theta_initial_guess=theta_start-L/b;
options=optimoptions('fsolve','Display','off');
theta=fsolve(arc_length_diff,theta_initial_guess, options);
end

% 定积分求解弧路径长度
function L = get_arc_length(theta_start,theta)
P=55;%螺距55cm
a=0;
b=P/(2*pi);%间距系数
arc_length = @(theta) sqrt(b^2+(a+b*theta).^2);
L=integral(arc_length,theta,theta_start);
end
