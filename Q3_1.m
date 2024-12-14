clc, clear, close all;
P=50;%从50开始往里面搜索
label=zeros(1,100);
for a=1:100
    t=200;
    X=zeros(1,224);
    Y=zeros(1,224);
    flag=0;
    while(flag==0)
        %使用包络线结果作为初始值
        [L,A,X,Y,V]=Q1_f(P,t);
        X=X';
        Y=Y';
        V=V';
        rectangles=[];
        for i=1:length(X)-1
            dx=X(i+1)-X(i);
            dy=Y(i+1)-Y(i);
            length_vec=sqrt(dx^2+dy^2);
            unit_vec=[dx, dy]/length_vec;
            perp_vec=[-unit_vec(2),unit_vec(1)];
            p1=[X(i),Y(i)]+perp_vec*15-unit_vec*27.5;
            p2=[X(i+1),Y(i+1)]+perp_vec*15+unit_vec*27.5;
            p3=[X(i+1),Y(i+1)]-perp_vec*15+unit_vec*27.5;
            p4=[X(i),Y(i)]-perp_vec*15-unit_vec*27.5;
            rectangles=[rectangles;struct('points',[p1; p2; p3; p4],'index',i)];
        end
        % 碰撞检测
        for i=1:length(rectangles)
            for j=i+2:length(rectangles)  % 检查非相邻的矩形
                if checkOverlap(rectangles(i).points,rectangles(j).points)
                    flag=1;
                end
            end
        end
        %如果没有碰撞 龙头还进入了调头区域，根据运动的连续性可得已经满足了 退出循环
        if(X(1)^2+Y(1)^2<=450^2)
            label(a)=t;
            break;
        end
        t=t+0.5;
    end
end
%% 极限圆进行碰撞检测
function [L, A, X, Y, V] = Q1_f(P,t)
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
X=R.*cos(theta);
Y=R.*sin(theta);
%包络线碰撞检测
    for i=2:5
        L(i-1)=sqrt((X(1)-X(i))^2+(Y(1)-Y(i))^2);
    end
lp_head=286;%板子参数
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
    if (beta1+beta2>=Angle)
        A=true;
    else
        A=false;
    end
end

%重叠检测
function overlap=checkOverlap(points1, points2)
%分离轴定理（SAT）来判断重叠
overlap=false;
%所有边作分离轴
edges=[points1(2,:)-points1(1,:);points1(3,:)-points1(2,:);points2(2,:)-points2(1,:);points2(3,:)-points2(2,:)];
for i=1:size(edges, 1)
    axis=[-edges(i,2),edges(i,1)];
    axis=axis/norm(axis);
    proj1=points1*axis';
    proj2=points2*axis';
    if max(proj1)<min(proj2)||max(proj2)<min(proj1)
        return;
    end
end
overlap=true; 
end
% 计算尾部节点角度的函数
function theta2=solve_theta2(theta1, X)
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