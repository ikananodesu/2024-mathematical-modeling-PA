clc, clear, close all;
% 加载路径数据
V=zeros(201,224);
i=0;
for t=-100:1:100
    i=i+1;
    V(i,:)=func(t);
end
V=V';
V=V/100;
surf(V);
shading interp

function V=func(t)
way=load("问题四路径.mat");
way=way.way; % 'way' 包含路径的 X 和 Y 坐标
% 出发点坐标
target_X=-396.085526860245;%初始X坐标
target_Y=213.579623119440;%初始Y坐标
r=241.1515;%圆半径
head=286;
body=165;
body_num=223;%板凳总量
points=zeros(2,body_num+1);%把手数量
label=zeros(1,body_num+1);%路径映射标签
len_points=zeros(1,body_num+1);%距离出发点的路径长
num_points=size(way, 2);%地图标签路径区间
V=zeros(1,body_num+1);%速度
a=zeros(1,body_num+1);%角速度
V(1)=100;%龙头速度
S=V(1)*t;%线路程
dist_to_target=sqrt((way(1,:)-target_X).^2+(way(2,:)-target_Y).^2);
[~,begin_index]=min(dist_to_target);%最近的路径标签
%% 求取龙头所在坐标
[lengthA,i]=func_get_len(S,begin_index,way);%出发点-龙头位置的弧积分等于运行距离
points(:,1)=way(:,i);%存储龙头坐标
len_points(1)=lengthA;
label(1)=i;
%% 回溯求解身体坐标
for i=1:body_num
    if(i==1)
        L=head;%龙头板
    else
        L=body;%其他板
    end
    dir=-1;%回溯
    [points(1,i+1),points(2,i+1),label(i+1)]=select_f(points(1,i),points(2,i),label(i),L,dir,way);
end
%% 分区规划计算速度
%10000 3988 10000

% 初始化极角 theta，用于极坐标表示下的速度更新
theta = zeros(1, body_num + 1);

% 计算每个点相对于 (0,0) 的极角 theta
for i = 1:body_num+1
    theta(i) = atan2d(points(2, i), points(1, i)); % 计算第 i 个点的极角，atan2d 返回角度，范围为 -180 到 180 度
end
d_theta=diff(theta);
A=find(d_theta<0);
for i=1:length(A)
    if(A(i)==223)
        d_theta(A(i))=d_theta(A(i)-1);
    elseif A(i)==1
        d_theta(A(i))=d_theta(A(i)+1);
    else
        d_theta(A(i))=(d_theta(A(i)+1)+d_theta(A(i)-1))/2;
    end
end
%% 更新速度迭代公式，基于极坐标
for i = 1:body_num
    if (label(i) <= 10000) % 入旋线
        V(i+1) = V(i) .* abs(cos((d_theta(i)) ./ 180)); % 根据极角的变化更新速度
    elseif (label(i) <= 13986 && label(i) > 10000) % 圆弧
        V(i+1) = V(i); % 保持速度不变
    else % 出旋线
        V(i+1) = V(i) .*abs(cos((d_theta(i)) ./ 180)); % 出旋线同样根据极角更新速度
    end
end

%% 绘制路径和长方形板凳
% figure;
% hold on;
% theta_circle = linspace(0,2*pi,2000);
% r_circle=450;%大圆半径
% x_circleA=r_circle * cos(theta_circle);
% y_circleA=r_circle * sin(theta_circle);
% plot(x_circleA, y_circleA,'m-','LineWidth', 1.5);
% plot(way(1,1:10000),way(2,1:10000),'r', 'LineWidth', 1.5)
% plot(way(1,10001:13988),way(2,10001:13988),'c-', 'LineWidth', 1.5)
% plot(way(1,13989:end),way(2,13989:end),'g', 'LineWidth', 1.5)
% scatter(points(1,2:end),points(2,2:end),20,'filled','o','MarkerFaceColor','b','MarkerEdgeColor','b'); % 绘制身体坐标点
% plot(points(1,:),points(2,:),'b.-'); % 绘制身体坐标点
% axis equal; % 确保x轴和y轴的比例相同
% grid on;

% 绘制相邻点之间的矩形
rectangles=[];
for i=1:length(points)-1
    dx=points(1, i+1)-points(1, i);
    dy=points(2, i+1)-points(2, i);
    length_vec=sqrt(dx^2+dy^2);
    unit_vec=[dx, dy]/length_vec;
    perp_vec=[-unit_vec(2), unit_vec(1)];
    p1=[points(1,i), points(2, i)]+perp_vec*15-unit_vec*27.5; % 上边
    p2=[points(1,i+1), points(2, i+1)]+perp_vec*15+unit_vec*27.5; % 上边
    p3=[points(1,i+1), points(2, i+1)]-perp_vec*15+unit_vec*27.5; % 下边
    p4=[points(1,i), points(2, i)]-perp_vec*15-unit_vec*27.5; % 下边
    rectangles=[rectangles;struct('points',[p1;p2;p3;p4],'index',i)];
    % fill([p1(1),p2(1),p3(1),p4(1)],[p1(2), p2(2), p3(2), p4(2)],'r','FaceAlpha',0.2,'EdgeColor','b');
end
% 碰撞检测：检查非相邻矩形是否重叠
for i=1:length(rectangles)
    for j=i+2:length(rectangles) % 检查非相邻的矩形
        if checkOverlap(rectangles(i).points,rectangles(j).points)
            % disp(['矩形',num2str(rectangles(i).index),'和',num2str(rectangles(j).index),'重叠']);
        end
    end
end
% 
% scatter(points(1,1),points(2,1),20,'filled','o','MarkerFaceColor','r','MarkerEdgeColor','r'); % 绘制身体坐标点
% xlabel('X坐标/cm');
% ylabel('Y坐标/cm');
% hold off;
end
%% 函数定义部分
% 输入量 X Y坐标 点的标签 板长 搜索方向 搜索标签点
function [X,Y,label_out] = select_f(points_x,points_y,label_in,L,dir,way)
    i=label_in;
    if (dir==1) %向前搜索
        while(i<length(way)&&abs(sqrt((way(1,i)-points_x)^2+(way(2,i)-points_y)^2)-L)>10)
            label_out=i;
            i=i+1;
        end
    else %向后搜索
        while(i>0&&abs(sqrt((way(1,i)-points_x)^2 + (way(2,i)-points_y)^2)-L)>10)
            label_out=i;
            i=i-1;
        end
    end
    X=way(1,label_out);
    Y=way(2,label_out);
end

%曲线积分反解路径点
function [length,i]=func_get_len(S,begin_index,way)
    %初始化
    num_points=size(way, 2);
    length=0;
    i=begin_index;% 出发点坐标
    if(S>0)%积分大于0 往前找
        while i<num_points&&length<S
            dX=way(1,i+1)-way(1,i);
            dY=way(2,i+1)-way(2, i);
            length=length+sqrt(dX^2+dY^2);
            i=i+1;
        end
    elseif (S<0)%积分小于0 往后找
        while i>1&&length<abs(S)
            dX=way(1,i-1)-way(1,i);
            dY=way(2,i-1)-way(2,i);
            length=length+sqrt(dX^2+dY^2);
            i=i-1;
        end
    end
end
%SAT算法检测碰撞
function overlap = checkOverlap(points1,points2)
    overlap=false;
    edges=[points1(2,:)-points1(1,:);points1(3,:)-points1(2,:);points2(2,:)-points2(1,:);points2(3,:)-points2(2,:)];
    for i=1:size(edges,1)
        axis=[-edges(i,2), edges(i,1)];
        axis=axis/norm(axis);
        proj1=points1*axis';
        proj2=points2*axis';
        if max(proj1)<min(proj2)||max(proj2)<min(proj1)
            return;
        end
    end
    overlap = true;
end