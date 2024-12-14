clc,clear,close;
% 加载路径数据
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
% 初始化时间范围 t 从 -100 到 100，间隔为 1
time_steps = -100:1:100;

% 创建一个矩阵用于存储每个时间步的速度
% 行数为 time_steps 的长度，列数为 body_num+1（包含龙头）
V_matrix = zeros(length(time_steps), body_num + 1);

% 对于每一个时间步长，更新龙头的位置和身体坐标，并记录速度
for t_idx = 1:length(time_steps)
    t = time_steps(t_idx);
    
    % 计算路径长度S, 假设速度 V(1) 是常量
    S = V(1) * t; % 这里你可以根据实际需要调整速度 V 和时间的关系

    % 求龙头的路径长度并确定坐标
    [lengthA, begin_index] = func_get_len(S, begin_index, way); % 使用更新的时间 t 来求解弧长
    points(:, 1) = way(:, begin_index); % 存储龙头新的坐标
    len_points(1) = lengthA;
    label(1) = begin_index;

    % 回溯求解身体的所有板凳位置
    for i = 1:body_num
        if(i == 1)
            L = head; % 龙头板长度
        else
            L = body; % 其他板的长度
        end
        dir = -1; % 回溯方向
        [points(1, i+1), points(2, i+1), label(i+1)] = select_f(points(1, i), points(2, i), label(i), L, dir, way);
    end

    % 计算当前时间步对应的每一块板的速度
    for i = 1:body_num
        if (label(i) <= 10000) % 入旋线
            V(i+1) = sqrt(points(1, i)^2 + points(2, i)^2) / sqrt(points(1, i+1)^2 + points(2, i+1)^2) * V(i);
        elseif (label(i) <= 13986 && label(i) > 10000) % 圆弧
            V(i+1) = V(i);
        else % 出旋线
            V(i+1) = sqrt(points(1, i)^2 + points(2, i)^2) / sqrt(points(1, i+1)^2 + points(2, i+1)^2) * V(i);
        end
    end

    % 将当前时间步的速度存储到矩阵中
    V_matrix(t_idx, :) = V;

end

% 最终的速度矩阵 V_matrix 包含每个时间步长下每一块板的速度
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