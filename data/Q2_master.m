% -------------------------------------------------------------------------
% 2020 Huawei Cup Mathematical Modeling: F
% Question 2 - master
% Author: Haoyu Miao
% Date: 2020/09/17
% Version: 2 Min-version
% This program is based on Yalmip + Gurobi solver
% -------------------------------------------------------------------------

% 从附件3中读取计划耗油速度数据 h(t) 以及 理想质心位置数据 c2(t)
clear
clc
h = xlsread('../data/附件3-问题2数据.xlsx',1,'B:B');
c2 = xlsread('../data/附件3-问题2数据.xlsx',2,'B:D');
h_split = zeros(120,1);
c2_split = zeros(120,3);
for i = 1:120
   h_split(i) = sum(h(60*(i-1)+1:60*i));
   c2_split(i,:) = c2(60*i,:); 
end


% -------------------------------------------------------------------------
% 常数 
% -- 油箱数量 
N = 6;
% -- 时间/s
T = 7200;
% -- 时间/min
T_split = 120;
% -- 供油速度上限/s
UB=[1.1,1.8,1.7,1.5,1.6,1.1];
% -- 大 M 法
M1=1000;
% 油箱长度
a = [1.5,2.2,2.4,1.7,2.4,2.4]';
% 油箱宽度
b = [0.9,0.8,1.1,1.3,1.2,1]';
% 油箱高度
c = [0.3,1.1,0.9,1.2,1,0.5]';
% 初始油量 （列向量）
V0 = [0.3,1.5,2.1,1.9,2.6,0.8]';



% 燃料密度
rho = 850;
% 飞行器质量
M = 3000;
% 油箱中心位置
P = [8.91304348,1.20652174,0.61669004;
     6.91304348,-1.39347826,0.21669004;
     -1.68695652,1.20652174,-0.28330996;
     3.11304348,0.60652174,-0.18330996;
     -5.28695652,-0.29347826,0.41669004;
     -2.08695652,-1.49347826,0.21669004;
];
% 飞行器（不载油）质心
c0 = [0,0,0];

% 初始油量的质量 （列向量）
m0 = V0 * rho;
% 最初各油箱的质量比例
m0_component = m0/sum(m0);


% 决策变量
% -- 1 是否供油
x = binvar(N,T_split);
% -- 2 供油速度,T_split意味着每分钟的，kg/min
y = sdpvar(N,T_split);
% -- 3 油的质量
m = sdpvar(N,T_split);
% -- 4 瞬时质心
co = sdpvar(N,T_split,3);
c1 = sdpvar(T_split,3);



% -------------------------------------------------------------------------
% 约束：
% -- 1 供油速度上限约束 与 下限约束
Constraints = [y>=0];
for i = 1:N
    % 每秒的
    % Constraints = [Constraints, y(i,:)<=UB(i)];
    % 每分的
    Constraints = [Constraints, y(i,:)<= 60 * UB(i)];
    fprintf('Constraints 1 - %d\n',i);
end

Constraints = [Constraints,m >= 0];
for t = 1:T_split
    Constraints = [Constraints,m(:,t) <= a.*b.*c*rho];
end
%%
% -- 2 供油持续时间约束
% for i = 1:N
%     for t = 1:T-60
%         Constraints = [Constraints, -M1*(1-(x(i,t+1))+x(i,t))<=sum(x(i,t+1:t+60))-60 <= M1*(1-(x(i,t+1))+x(i,t))];
%     end
%     fprintf('Constraints 2 - %d\n',i);
% end

% -- 3 油箱数约束
for t = 1:T_split
    Constraints = [Constraints,sum(x(2:5,t))<=2];
    Constraints = [Constraints,sum(x(1:6,t))<=3];
end

%%
% -- 4 耗油量约束
for i = 1:N
    Constraints = [Constraints,-M1*x(i,:)<=y(i,:)<=M1*x(i,:)]; 
    fprintf('Constraints 4 - %d\n',i);
end
Constraints = [Constraints,sum(y(2:5,:))>=h_split'];

%%
% -- 5 油的初始质量约束
Constraints = [Constraints,m(:,1) - m0==0];

%%
% -- 6 油的动态方程约束
for t = 1:T_split-1
    Constraints = [Constraints,m(1,t+1)==m(1,t) - y(1,t)];
    Constraints = [Constraints,m(2,t+1)==m(2,t) - y(2,t) + y(1,t)];
    Constraints = [Constraints,m(3,t+1)==m(3,t) - y(3,t)];
    Constraints = [Constraints,m(4,t+1)==m(4,t) - y(4,t)];
    Constraints = [Constraints,m(5,t+1)==m(5,t) - y(5,t) + y(6,t)];
    Constraints = [Constraints,m(6,t+1)==m(6,t) - y(6,t)];
    fprintf('Constraints 6 - %d\n',t);
end

%%
% -- 7 油体瞬时质心约束
% -- 如果有俯仰角，需要改变这组约束
for t = 1:T_split
    Constraints = [Constraints,co(:,t,1) == P(:,1)]; 
    Constraints = [Constraints,co(:,t,2) == P(:,2)]; 
    for i = 1:N
        Constraints = [Constraints,co(i,t,3) == P(i,3) - c(i)/2 + m(i,t)/(2*a(i)*b(i)*rho)]; 
    end 
    fprintf('Constraints 7 - %d\n',t);
end

%% 
% -- 8 飞行器瞬时质心动态方程
% 这里我需要一个估计质量
M_estimator = zeros(T_split,1);
M_estimator(1) = M + sum(m0);
for t = 2:T_split
   M_estimator(t) = M_estimator(t-1) - h_split(t-1);
end

m_estimator = zeros(N,T_split);
m_estimator(:,1) = m0;
for t = 2:T_split
   m_estimator(:,t) = m_estimator(:,t-1) - h_split(t) * m0_component;
end
constant1 = zeros(N,T_split);
for i = 1:N
    constant1(i,:) = m_estimator(i,:)/(2*a(i)*b(i)*rho);
end

%%
for t = 1:T_split
    % -- x
    Constraints = [Constraints,M_estimator(t) * c1(t,1)== m(:,t)'* P(:,1)]; 
    % -- y
    Constraints = [Constraints,M_estimator(t) * c1(t,2)== (m(:,t)'* P(:,2) + M * c0(2))]; 
    % -- z
    Constraints = [Constraints,M_estimator(t) * c1(t,3)== m(:,t)'* (P(:,3) - c/2) + m(:,t)'* constant1(:,t) + M * c0(3)];
    fprintf('Constraints 8 - %d\n',t);
end


fprintf('All constraints added\n');
% 目标函数
% #########################################################################
% 计算 c1 的表达式
error = c1-c2_split;
Objective = max(sqrt(sum((error).*(error),2)));
% Objective = -sum(sum(x));
% #########################################################################

% 求解
options = sdpsettings('solver','gurobi');
sol = optimize(Constraints,Objective,options);

% c1_val = value(c1);
% plot(c1_val(:,1));
% hold on
% plot(c2_split(:,1));





























