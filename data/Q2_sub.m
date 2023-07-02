% -------------------------------------------------------------------------
% 2020 Huawei Cup Mathematical Modeling: F
% Question 2 - sub
% Author: Haoyu Miao
% Date: 2020/09/20
% Version: 2 Min-version
% This program is based on Yalmip + Gurobi solver
% -------------------------------------------------------------------------

% x 是 6 * 120 的 binary变量
% 对于第 i 分钟
Q2_master;
%%
c2_val = c2;
Obj_val = value(Objective);
co_val = value(co);
m_val = value(m);
x_val = value(x);
y_val = value(y);
c1_val = value(c1);

% -------------------------------------------------------------------------
% 常数 
% -- 油箱数量 
N = 6;
% -- 时间/s
% T = 7200;
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

m_init = m0;
y_all = [];
m_all = [];
c1_all = [];

A = [-1,0,0,0,0,0;
         1,-1,0,0,0,0;
         0,0,-1,0,0,0;
         0,0,0,-1,0,0;
         0,0,0,0,-1,1;
         0,0,0,0,0,-1];

%%
for ii = 1:120    
	T=60;
	y = sdpvar(6,60);
	m = sdpvar(6,60);
    co = sdpvar(N,T,3);
    c1 = sdpvar(T,3);
	Constraints = [];

	% 上下限约束
	Constraints = [Constraints,y>=0];
	for j = 1:6
		Constraints = [Constraints,y(j,:)<=UB(j)];	
    end

%     Constraints = [Constraints,m >= 0];
%     for t = 1:T
%         Constraints = [Constraints,m(:,t) <= a.*b.*c*rho];
%     end
    
	% 总供油量约束
	for j = 1:6
		if y_val(j,ii) == 0
			Constraints = [Constraints, y(j,:)==0];	
        else
			Constraints = [Constraints, sum(y(j,:))<=y_val(j,ii)];
        end
    end

    Constraints = [Constraints,h(60*(ii-1)+1:60*ii)'<=sum(y(2:5,:))];
	% 质量初始约束，m_init 要更新
    Constraints = [Constraints, m(:,1)==m_init]; 
	
	% 油的动态方程约束
	for t = 1:T-1
    	Constraints = [Constraints,m(1,t+1)==m(1,t) - y(1,t)];
 		Constraints = [Constraints,m(2,t+1)==m(2,t) - y(2,t) + y(1,t)];
    	Constraints = [Constraints,m(3,t+1)==m(3,t) - y(3,t)];
    	Constraints = [Constraints,m(4,t+1)==m(4,t) - y(4,t)];
   		Constraints = [Constraints,m(5,t+1)==m(5,t) - y(5,t) + y(6,t)];
    	Constraints = [Constraints,m(6,t+1)==m(6,t) - y(6,t)];
	end
	
	% -- 7 油体瞬时质心约束
	% -- 如果有俯仰角，需要改变这组约束
	for t = 1:T
    		Constraints = [Constraints,co(:,t,1) == P(:,1)]; 
    		Constraints = [Constraints,co(:,t,2) == P(:,2)]; 
    		for i = 1:N
        		Constraints = [Constraints,co(i,t,3) == P(i,3) - c(i)/2 + m(i,t)/(2*a(i)*b(i)*rho)]; 
    		end 
    		fprintf('Constraints 7 - %d\n',t);
    end

    
% 质量估计
    % 估计总质量
    M_estimator = zeros(T,1);
    M_estimator(1) = M + sum(m_init);
    for t = 2:T
        M_estimator(t) = M_estimator(t-1) - sum(y_val(2:5,ii))/60;
    end    
    
    m_estimator = zeros(N,T);
    m_estimator(:,1) = m_init;
    
    for t = 2:T
        m_estimator(:,t) = m_estimator(:,t-1) + A * y_val(:,ii)/60;
    end
    
    constant1 = zeros(N,T);
    for i = 1:N
        constant1(i,:) = m_estimator(i,:)/(2*a(i)*b(i)*rho);
    end

    for t = 1:T
        % -- x
        Constraints = [Constraints,M_estimator(t) * c1(t,1)== m(:,t)'* P(:,1)]; 
        % -- y
        Constraints = [Constraints,M_estimator(t) * c1(t,2)== m(:,t)'* P(:,2)]; 
        % -- z
        Constraints = [Constraints,M_estimator(t) * c1(t,3)== m(:,t)'* (P(:,3) - c/2) + m(:,t)'* constant1(:,t)];
        fprintf('Constraints 8 - %d\n',t);
    end

    fprintf('All constraints added');
    % 目标函数
    % #########################################################################
    % 计算 c1 的表达式
    error = c1-c2((ii-1)*60+1:ii*60,:);
    Objective = max(sqrt(sum((error).*(error),2)));
    % Objective = -sum(sum(x));
    % #########################################################################

    % 求解
    options = sdpsettings('solver','gurobi');
    sol = optimize(Constraints,Objective,options);
    
    y_all = [y_all,value(y)];
    m_all = [m_all,value(m)];
    % co_all = [co_all;value(co)];
    c1_all = [c1_all;value(c1)];
    fprintf('正在计算%d 分钟\n',ii);
    
    % 更新部分
    m_init = m_val(:,ii+1);
end