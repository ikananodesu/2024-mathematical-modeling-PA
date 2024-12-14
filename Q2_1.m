clc,clear,close all;
%% 包络线方法求解

%搜索器
% t=441;
% %极限距离：
% Lmin=2*sqrt((27.5+5.5/2)^2+15^2);%构建的最小包络线
% for t=410:0.01:430
%     [L,A,X,Y]=Q1_f(t);
%     if(min(L(2:end))<Lmin)
%         break;
%     end
%     if(A)
%         break;
%     end
% end
t=413.16;
[L,A,X,Y,V]=Q1_f(t);
%转单位
X=X'/100;
Y=Y'/100;
V=V'/100;

%% 极限圆进行碰撞检测
function [L,A,X,Y,V]=Q1_f(t)

%定义一下螺线，以螺心作为0°点，每逆时针旋转360°代表往外转一圈
%定义阿基米德螺线情况
%r=a+b*sita
P=55;%螺距55cm
a=0;
b=P/(2*pi);%间距系数
body_num=223;%板凳总量
theta=zeros(1,body_num+1);%初始化极坐标存储
V=100;%龙头速度100cm/s
S=V*t;%龙头移动弧长
theta_start=16*2*pi;%第16圈出发
theta(1)=solve_theta_from_arc_length(S,theta_start);%初始化龙头位置
%% 定义龙头和身体的属性
% 板长 用来定义头节点i1和龙头尾节点i2的距离L
lp_head=286;
lp_body=165;
%1个龙头221个龙身1个龙尾 一共224个把手需要推导
for i=1:body_num
    if i==1
        theta(i+1)=solve_theta2(theta(i),lp_head);
    else
        theta(i+1)=solve_theta2(theta(i),lp_body);
    end
end
% 每个点对应于中心的距离
R=get_length(theta);
V=zeros(1,body_num+1);
V(1)=100;
for i=1:body_num
    V(i+1)=V(i)*cos((theta(i+1)-theta(i))/180);
end
X=R.*cos(theta);
Y=R.*sin(theta);
figure(1);
hold on;
scatter(X,Y);
theta_values=linspace(0,theta_start,1000);
R_values=get_length(theta_values);
X_values=R_values.*cos(theta_values);
Y_values=R_values.*sin(theta_values);
plot(X_values,Y_values,'r','LineWidth',1.5);
plot(X,Y,'b.-')
%% 判断是否碰撞
for i=2:5
    L(i-1)=sqrt((X(1)-X(i))^2+(Y(1)-Y(i))^2);
end
lp_head=286;
lp_body=165;
beta1=atand(15/lp_head);
beta2=atand(15/lp_body);
AB=[X(1)-X(2),Y(1)-Y(2)];
BC=[X(2)-X(3),Y(2)-Y(3)];
A=dot(AB,BC);
mAB=norm(AB);
mBC=norm(BC);
Angle=A/(mAB*mBC);
Angle=acosd(Angle);
if(beta1+beta2>=Angle)
    A=true;
else
    A=false;
end
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
function theta=solve_theta_from_arc_length(L,theta_start)
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