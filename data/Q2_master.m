% -------------------------------------------------------------------------
% 2020 Huawei Cup Mathematical Modeling: F
% Question 2 - master
% Author: Haoyu Miao
% Date: 2020/09/17
% Version: 2 Min-version
% This program is based on Yalmip + Gurobi solver
% -------------------------------------------------------------------------

% �Ӹ���3�ж�ȡ�ƻ������ٶ����� h(t) �Լ� ��������λ������ c2(t)
clear
clc
h = xlsread('../data/����3-����2����.xlsx',1,'B:B');
c2 = xlsread('../data/����3-����2����.xlsx',2,'B:D');
h_split = zeros(120,1);
c2_split = zeros(120,3);
for i = 1:120
   h_split(i) = sum(h(60*(i-1)+1:60*i));
   c2_split(i,:) = c2(60*i,:); 
end


% -------------------------------------------------------------------------
% ���� 
% -- �������� 
N = 6;
% -- ʱ��/s
T = 7200;
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


% ���߱���
% -- 1 �Ƿ���
x = binvar(N,T_split);
% -- 2 �����ٶ�,T_split��ζ��ÿ���ӵģ�kg/min
y = sdpvar(N,T_split);
% -- 3 �͵�����
m = sdpvar(N,T_split);
% -- 4 ˲ʱ����
co = sdpvar(N,T_split,3);
c1 = sdpvar(T_split,3);



% -------------------------------------------------------------------------
% Լ����
% -- 1 �����ٶ�����Լ�� �� ����Լ��
Constraints = [y>=0];
for i = 1:N
    % ÿ���
    % Constraints = [Constraints, y(i,:)<=UB(i)];
    % ÿ�ֵ�
    Constraints = [Constraints, y(i,:)<= 60 * UB(i)];
    fprintf('Constraints 1 - %d\n',i);
end

Constraints = [Constraints,m >= 0];
for t = 1:T_split
    Constraints = [Constraints,m(:,t) <= a.*b.*c*rho];
end
%%
% -- 2 ���ͳ���ʱ��Լ��
% for i = 1:N
%     for t = 1:T-60
%         Constraints = [Constraints, -M1*(1-(x(i,t+1))+x(i,t))<=sum(x(i,t+1:t+60))-60 <= M1*(1-(x(i,t+1))+x(i,t))];
%     end
%     fprintf('Constraints 2 - %d\n',i);
% end

% -- 3 ������Լ��
for t = 1:T_split
    Constraints = [Constraints,sum(x(2:5,t))<=2];
    Constraints = [Constraints,sum(x(1:6,t))<=3];
end

%%
% -- 4 ������Լ��
for i = 1:N
    Constraints = [Constraints,-M1*x(i,:)<=y(i,:)<=M1*x(i,:)]; 
    fprintf('Constraints 4 - %d\n',i);
end
Constraints = [Constraints,sum(y(2:5,:))>=h_split'];

%%
% -- 5 �͵ĳ�ʼ����Լ��
Constraints = [Constraints,m(:,1) - m0==0];

%%
% -- 6 �͵Ķ�̬����Լ��
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
% -- 7 ����˲ʱ����Լ��
% -- ����и����ǣ���Ҫ�ı�����Լ��
for t = 1:T_split
    Constraints = [Constraints,co(:,t,1) == P(:,1)]; 
    Constraints = [Constraints,co(:,t,2) == P(:,2)]; 
    for i = 1:N
        Constraints = [Constraints,co(i,t,3) == P(i,3) - c(i)/2 + m(i,t)/(2*a(i)*b(i)*rho)]; 
    end 
    fprintf('Constraints 7 - %d\n',t);
end

%% 
% -- 8 ������˲ʱ���Ķ�̬����
% ��������Ҫһ����������
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
% Ŀ�꺯��
% #########################################################################
% ���� c1 �ı��ʽ
error = c1-c2_split;
Objective = max(sqrt(sum((error).*(error),2)));
% Objective = -sum(sum(x));
% #########################################################################

% ���
options = sdpsettings('solver','gurobi');
sol = optimize(Constraints,Objective,options);

% c1_val = value(c1);
% plot(c1_val(:,1));
% hold on
% plot(c2_split(:,1));





























