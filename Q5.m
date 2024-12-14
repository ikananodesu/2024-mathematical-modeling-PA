clc, clear, close all;
%搜索每个速度于-100到100时间窗内的各峰值速度
time=-200:1:300;
Vbegin=1.99:0.001:2;
Vselect={length(time),length(Vbegin)};
Aselect={length(time),length(Vbegin)};
label=zeros(length(time),length(Vbegin));
for i=1:length(time)
    for j=1:length(Vbegin)
        [Vselect{i,j},Aselect{i,j}]=get_vmax(Vbegin(j),time(i));
        if(max(Vselect{i,j})>=2)
            label(i,j)=1;
        else
            label(i,j)=0;
        end
    end
end
% 
% %搜索每个速度于-100到100时间窗内的各峰值速度
% time=-100:1:100;
% Vbegin=1;
% Vselect=zeros(length(time),224);
% Aselect=zeros(length(time),224);
% label=zeros(length(time),length(Vbegin));
% for i=1:length(time)
%     for j=1:length(Vbegin)
%         [Vselect(i,:),Aselect(i,:)]=get_vmax(Vbegin(j),time(i));
%     end
% end
% Vselect=Vselect';
function [V,a]=get_vmax(Vbegin,time)
% 加载路径数据
way=load("问题四路径.mat");
way=way.way;
target_X=-396.085526860245;%初始X坐标
target_Y=213.579623119440;%初始Y坐标
r=241.1515;%圆半径
head=286;
body=165;
body_num=223;%板凳总量
points=zeros(2,body_num+1);%把手
label=zeros(1,body_num+1);%映射路径标签i
len_points=zeros(1,body_num+1);
num_points=size(way,2);%搜索点数量
V=zeros(1,body_num+1);
a=zeros(1,body_num+1);
V(1)=Vbegin;%龙头速度
t=time;%运行时间
S=V(1)*t;%弧线路程
dist_to_target=sqrt((way(1,:)-target_X).^2+(way(2,:)-target_Y).^2);
[~,begin_index]=min(dist_to_target);%映射标签
%% 求取龙头所在坐标
% 曲线积分求解得到目标点所在坐标
[lengthA,i]=func_get_len(S,begin_index,way);
% 存储龙头位置的坐标
points(:,1)=way(:,i);
len_points(1)=lengthA;
label(1)=i;
%% 回溯求解身体坐标
for i=1:body_num
    if(i==1)
        L=head;
    else
        L=body;
    end
    dir=-1;%反方向搜索
    [points(1,i+1), points(2,i+1),label(i+1)]=select_f(points(1,i),points(2,i),label(i),L,dir,way);
end

% 初始化极角 theta，用于极坐标表示下的速度更新
theta = zeros(1, body_num + 1);

% 计算每个点相对于 (0,0) 的极角 theta
for i = 1:body_num+1
    theta(i) = atan2d(points(2, i), points(1, i)); % 计算第 i 个点的极角，atan2d 返回角度，范围为 -180 到 180 度
end
d_theta=diff(theta);
A=find(d_theta<0);
for ii=1:length(A)
    if(A(ii)==223)
        d_theta(A(ii))=d_theta(A(ii)-1);
    elseif A(ii)==1
        d_theta(A(ii))=d_theta(A(ii)+1);
    else
        d_theta(A(ii))=(d_theta(A(ii)+1)+d_theta(A(ii)-1))/2;
    end
end
%% 更新速度迭代公式，基于极坐标
for i = 1:body_num
    if (label(i) <= 10000) % 入旋线
        V(i+1) = V(i) .* abs(cos((d_theta(i)) ./ 180)); % 根据极角的变化更新速度
    elseif (label(i) <= 13986 && label(i) > 10000) % 圆弧
        V(i+1) = V(i); % 保持速度不变
    else % 出旋线
        V(i+1) = V(i) ./(abs(cos((d_theta(i)) ./ 180))); % 出旋线同样根据极角更新速度
    end
end
% 绘制路径和点
% figure;
% hold on;
% plot(way(1,:),way(2,:)); % 绘制路径
% scatter(points(1,:),points(2,:)); % 绘制身体坐标点
end
%% 函数定义部分
%输入量 X Y坐标 点的标签 板长 搜索方向 输出所在点的标签
function [X,Y,label_out]=select_f(points_x, points_y, label_in,L,dir,way)
%从出发点翻转搜索求解最近点
i=label_in;
if(dir==1)%向前搜索
    while(i<length(way)&&abs(sqrt((way(1,i)-points_x)^2+(way(2,i)-points_y)^2)-L)>10)
        label_out=i;
        i=i+1;
    end
else%向后搜索
    while(i>0&&abs(sqrt((way(1,i)-points_x)^2+(way(2,i)-points_y)^2)-L)>10)
        label_out=i;
        i=i-1;
    end
end
X=way(1,label_out);
Y=way(2,label_out);
end

% 定积分反解路径
function [length,i]=func_get_len(S,begin_index,way) % 定积分值和对应的路径编号
% 积分初始化
num_points=size(way,2);
length=0;
i=begin_index; % 出发点坐标
if(S>0) % 如果S大于0代表需要往前推导龙头所在位置
    while i<num_points&&length<S
        dX=way(1,i+1)-way(1, i);
        dY=way(2,i+1)-way(2, i);
        length=length+sqrt(dX^2+dY^2);
        i=i+1;
    end
elseif(S<0) % 如果S小于0代表需要往后推导龙头所在位置
    while i>1&&length<abs(S) % 使用 abs(S) 来处理负值
        dX=way(1,i-1)-way(1,i);
        dY=way(2,i-1)-way(2,i);
        length=length+sqrt(dX^2+dY^2);
        i=i-1;
    end
end
end