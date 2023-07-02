% -------------------------------------------------------------------------
% 2020 Huawei Cup Mathematical Modeling: F
% Question 2 - sub
% Author: Haoyu Miao
% Date: 2020/09/20
% Version: 2 Min-version
% This program is based on Yalmip + Gurobi solver
% -------------------------------------------------------------------------

% x �� 6 * 120 �� binary����
% ���ڵ� i ����
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
% ���� 
% -- �������� 
N = 6;
% -- ʱ��/s
% T = 7200;
% -- ʱ��/min
T_split = 120;
% -- �����ٶ�����/s
UB=[1.1,1.8,1.7,1.5,1.6,1.1];
% -- �� M ��
M1=1000;
% ���䳤��
a = [1.5,2.2,2.4,1.7,2.4,2.4]';
% ������
b = [0.9,0.8,1.1,1.3,1.2,1]';
% ����߶�
c = [0.3,1.1,0.9,1.2,1,0.5]';
% ��ʼ���� ����������
V0 = [0.3,1.5,2.1,1.9,2.6,0.8]';
% ȼ���ܶ�
rho = 850;
% ����������
M = 3000;
% ��������λ��
P = [8.91304348,1.20652174,0.61669004;
     6.91304348,-1.39347826,0.21669004;
     -1.68695652,1.20652174,-0.28330996;
     3.11304348,0.60652174,-0.18330996;
     -5.28695652,-0.29347826,0.41669004;
     -2.08695652,-1.49347826,0.21669004;
];
% �������������ͣ�����
c0 = [0,0,0];

% ��ʼ���������� ����������
m0 = V0 * rho;
% ������������������
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

	% ������Լ��
	Constraints = [Constraints,y>=0];
	for j = 1:6
		Constraints = [Constraints,y(j,:)<=UB(j)];	
    end

%     Constraints = [Constraints,m >= 0];
%     for t = 1:T
%         Constraints = [Constraints,m(:,t) <= a.*b.*c*rho];
%     end
    
	% �ܹ�����Լ��
	for j = 1:6
		if y_val(j,ii) == 0
			Constraints = [Constraints, y(j,:)==0];	
        else
			Constraints = [Constraints, sum(y(j,:))<=y_val(j,ii)];
        end
    end

    Constraints = [Constraints,h(60*(ii-1)+1:60*ii)'<=sum(y(2:5,:))];
	% ������ʼԼ����m_init Ҫ����
    Constraints = [Constraints, m(:,1)==m_init]; 
	
	% �͵Ķ�̬����Լ��
	for t = 1:T-1
    	Constraints = [Constraints,m(1,t+1)==m(1,t) - y(1,t)];
 		Constraints = [Constraints,m(2,t+1)==m(2,t) - y(2,t) + y(1,t)];
    	Constraints = [Constraints,m(3,t+1)==m(3,t) - y(3,t)];
    	Constraints = [Constraints,m(4,t+1)==m(4,t) - y(4,t)];
   		Constraints = [Constraints,m(5,t+1)==m(5,t) - y(5,t) + y(6,t)];
    	Constraints = [Constraints,m(6,t+1)==m(6,t) - y(6,t)];
	end
	
	% -- 7 ����˲ʱ����Լ��
	% -- ����и����ǣ���Ҫ�ı�����Լ��
	for t = 1:T
    		Constraints = [Constraints,co(:,t,1) == P(:,1)]; 
    		Constraints = [Constraints,co(:,t,2) == P(:,2)]; 
    		for i = 1:N
        		Constraints = [Constraints,co(i,t,3) == P(i,3) - c(i)/2 + m(i,t)/(2*a(i)*b(i)*rho)]; 
    		end 
    		fprintf('Constraints 7 - %d\n',t);
    end

    
% ��������
    % ����������
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
    % Ŀ�꺯��
    % #########################################################################
    % ���� c1 �ı��ʽ
    error = c1-c2((ii-1)*60+1:ii*60,:);
    Objective = max(sqrt(sum((error).*(error),2)));
    % Objective = -sum(sum(x));
    % #########################################################################

    % ���
    options = sdpsettings('solver','gurobi');
    sol = optimize(Constraints,Objective,options);
    
    y_all = [y_all,value(y)];
    m_all = [m_all,value(m)];
    % co_all = [co_all;value(co)];
    c1_all = [c1_all;value(c1)];
    fprintf('���ڼ���%d ����\n',ii);
    
    % ���²���
    m_init = m_val(:,ii+1);
end