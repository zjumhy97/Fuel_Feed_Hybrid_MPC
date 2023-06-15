# bb_tree_set_cover.py
# author: Haoyu Miao
# date: 2022/03/02

import gurobipy as gp
from gurobipy import GRB
from operator import attrgetter
import copy
from heapq import *
import numpy as np
import time
import itertools
counter = itertools.count() # create a iterator from 0

class BBTree():
    def __init__(self, frontier=set(), tolerance=1e-3, treeUB=float('inf'), treeLB=-float('inf')):
        # 按照 paper 中的说法，frontier，其中每个集合的下界，以及 treeUB，这三个是必须的
        # 【重要！！！】 frontier 中每个 node 的下界应该是 parent node 对应对偶问题的界
        self.frontier = frontier
        self.tolerance = tolerance
        self.treeUB = treeUB
        self.treeLB = treeLB
        self.bestnode = None

    def selectSubproblem(self, heap):
        # 从 heap 中选取一个 node，选取的规则是：best-first，哪个下界最小，选哪个
        # such that UB-LB > self.tolerance
        node = min(heap, key=attrgetter('LB'))
        return node

    def bbtreeSolve(self):
        # 堆 heap 是用来控制循环的机制，heap 不为空意味着还有可以进行分枝的节点，算法仍要继续
        # 但是 heap 为空意味着所有的节点都被剪枝剪掉了，不需要
        # 构建 heap 的时候应该把 frontier 中的节点都纳入进来
        heap = []
        for node in self.frontier:
            # heappush(heap, node)
            heap.append(node)
        # print(heap)
        # Q: 根节点需要求解吗？
        # Q: 需要更新B&B树上界和最新节点吗？
        nodecount = 0
        while len(heap) > 0: # 这个循环是边界 frontier 迭代的循环，第0次，frontier就是根节点
            print('-----------  ', nodecount, '  -----------')
            # 用 len(heap) 作为循环条件的缺点——每次求解bb树都要求解完，而不是中途求到feasible的解也可以退出
            nodecount += 1
            # print("Heap Size: ", len(heap))
            # solve the subproblem chosen by selectSubproblem function, get relaxation value \theta(V^{i})
            # 满足条件的话，从边界中选一个节点求解, 2:optimal 3:infeasible 5:unbounded
            node = self.selectSubproblem(heap) # 从 frontier 中还是 heap 中？
            # 选择一个节点后，应该把节点从 heap 中弹出
            heap.remove(node) # heap 中移除，但是 frontier 中不移除

            # buildProblem 这里需要参数的哦~
            prob = node.buildProblem()
            prob.optimize()
            print('node:', node.v_underline, node.v_bar)
            if prob.Status == 2: # the status is optimal, not infeasible (3) or unbounded (5)
                res = prob.objVal
                if res >= self.treeUB - self.tolerance: # pruning
                    print("Relaxed Problem Stinks. Killing this branch.")
                    node.LB = res
                    pass
                elif node.isIntegral(prob): # solution update
                    print("New Best Integral solution. Optimal = " ,res)
                    self.treeUB = res
                    node.LB = res
                    node.UB = res
                    self.bestnode = node

                    # break 是为了测试最快的速度
                    # break
                else: # branching
                    # 现在开始处理分支，争取2022/03/05 23:30前处理完
                    newnodes = node.branch(model=prob)
                    self.frontier.remove(node)
                    for newnode in newnodes:
                        self.frontier.append(newnode)
                        heap.append(newnode)
                        # heappush(self.frontier, newnode)
                        # heappush(heap, newnode)
        # print("Nodes searched: ", nodecount)
        return self.treeUB, self.bestnode, self.frontier
        # 每次求解完需要返回哪些东西传给下一个B&B树？
        # 每次求解完后不需要的节点所占的内存如何释放？

class BBTreeNode():
    def __init__(self, vars=set(), bool_vars=None, objective=0, constraints=[], UB=float('inf'),
                 LB=-float('inf'),horizon_is_1=True, v_underline=[], v_bar=[]):
        self.vars = vars
        self.bool_vars = bool_vars
        self.objective = objective
        self.constraints = constraints
        self.UB = UB
        self.LB = LB
        self.horizon_is_1 = horizon_is_1
        # v_underline 和 v_bar 就按一维数组来实现
        self.v_underline = v_underline
        self.v_bar = v_bar
        self.children = [] # Q: what is the use of children?

    def buildProblem(self):
    # def buildProblem(self, x0, ):
        # 所以这里实际要求的是原问题 P 的松弛问题 P(V) 的对偶问题 D(V)，实验证明对偶问题求解成功
        # Q: buildProblem 需要什么参数？
        # 如果满足强对偶性的要求，求解得到的 D(V) 的最优值和 P(V)的最优值，也就是节点的下界 node.LB
        # prob = grb.Model('subproblem')
        # 添加对偶变量（对偶乘子），目标函数，约束
        # --------------------------------------------------------------------------------------------------------------
        # 以下内容是特定问题的
        # 参数需不需要移到class的外面？因为每次赋值都会占用一定的时间。

        T = 2  # time horizon
        u_dimension = 7  # input dimension
        u_binary_dimension = 4  # input dimension
        prob_x0 = np.array([0.3, 0.4, 0.5, 0.6])
        prob_Q = np.identity(4)  # Q^TQ=Q
        prob_QT = np.transpose(prob_Q)
        prob_R = np.zeros((7, 7))
        prob_R[0][0] = 1
        prob_RT = np.transpose(prob_R)
        prob_A = np.identity(4)
        prob_AT = np.transpose(prob_A)
        prob_B = np.array([[0.1, 0.2, 0.3, 1, 0, 0, 0],
                       [0.1, 0.2, 0.3, 0, 1, 0, 0],
                       [0.1, 0.2, 0.3, 0, 0, 1, 0],
                       [0.1, 0.2, 0.3, 0, 0, 0, 1]])
        prob_BT = np.transpose(prob_B)
        prob_F = np.identity(4)
        prob_FT = np.transpose(prob_F)
        prob_G = np.array([[0, 0, 0, 1, 0, 0, 0],
                       [0, 0, 0, 0, 1, 0, 0],
                       [0, 0, 0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 0, 0, 1]])
        prob_GT = np.transpose(prob_G)
        prob_V = np.array([[0, 0, 0, 1, 0, 0, 0],
                       [0, 0, 0, 0, 1, 0, 0],
                       [0, 0, 0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 0, 0, 1]])
        prob_VT = np.transpose(prob_V)
        # prob_v_bar = np.array([1, 1, 1, 1])
        # prob_v_underline = np.array([0, 0, 0, 0])
        prob_h = np.array([100, 100, 100, 100])
        horizon_is_1 = self.horizon_is_1
        if horizon_is_1: # T=1
            # 二元变量线性松弛的上下界
            # 在 B&B 框架中，最重要的应该是下面这两个变量，这两个变量是二元变量的上下界，这两个应该是传进来的
            # prob_v_bar = np.array([1, 1, 1, 1])
            # prob_v_underline = np.array([0, 0, 0, 0])
            prob_v_bar = self.v_bar
            prob_v_underline = self.v_underline

            model = gp.Model('miqp_dual')
            model.setParam('Method', 1)
            model.setParam('outPutFlag', 0)

            prob_lambda_d = model.addMVar((4, T + 1), vtype=GRB.CONTINUOUS, name='lambda_d', lb=-1e20, ub=1e20)
            prob_mu_d = model.addMVar((4, T), vtype=GRB.CONTINUOUS, name='mu_d', lb=0, ub=1000)
            prob_nu_underline_d = model.addMVar((4, T), vtype=GRB.CONTINUOUS, name='nu_underline_d', lb=0, ub=1e20)
            prob_nu_bar_d = model.addMVar((4, T), vtype=GRB.CONTINUOUS, name='nu_bar_d', lb=0, ub=1e20)
            prob_rou_d = model.addMVar((4, T + 1), vtype=GRB.CONTINUOUS, name='rou_d', lb=-1e20, ub=1e20)
            prob_delta_d = model.addMVar((7, T), vtype=GRB.CONTINUOUS, name='delta_d', lb=-1e20, ub=1e20)

            model.setObjective(-0.25 * sum(prob_rou_d[:, t] @ prob_rou_d[:, t] for t in range(T + 1))
                           - sum(0.25 * (prob_delta_d[:, t] @ prob_delta_d[:, t]) + prob_mu_d[:, t] @ prob_h
                                 + prob_nu_bar_d[:, t] @ prob_v_bar - prob_nu_underline_d[:, t] @ prob_v_underline
                                 for t in range(T)) - prob_lambda_d[:, 0] @ prob_x0, GRB.MAXIMIZE)

            for t in range(T):
                model.addConstr(prob_QT @ prob_rou_d[:, t] + prob_lambda_d[:, t] - prob_AT @ prob_lambda_d[:, t + 1]
                            + prob_FT @ prob_mu_d[:, t] == 0, name='d1')

            model.addConstr(prob_QT @ prob_rou_d[:, T] + prob_lambda_d[:, T] == 0, name='d2')

            for t in range(T):
                model.addConstr(prob_RT @ prob_delta_d[:, t] - prob_BT @ prob_lambda_d[:, t + 1]
                            + prob_GT @ prob_mu_d[:, t] + prob_VT @ prob_nu_bar_d[:, t]
                            - prob_VT @ prob_nu_underline_d[:, t] == 0, name='d3')
        else: # T > 1
            # 二元变量线性松弛的上下界
            # 在 B&B 框架中，最重要的应该是下面这两个变量，这两个变量是二元变量的上下界，这两个应该是传进来的

            prob_v_bar = np.array(self.v_bar).reshape(T, u_binary_dimension).T
            prob_v_underline = np.array(self.v_underline).reshape(T, u_binary_dimension).T

            model = gp.Model('miqp_dual')
            model.setParam('Method', 1)
            model.setParam('outPutFlag', 0)

            prob_lambda_d = model.addMVar((4, T + 1), vtype=GRB.CONTINUOUS, name='lambda_d', lb=-1e20, ub=1e20)
            prob_mu_d = model.addMVar((4, T), vtype=GRB.CONTINUOUS, name='mu_d', lb=0, ub=1000)
            prob_nu_underline_d = model.addMVar((4, T), vtype=GRB.CONTINUOUS, name='nu_underline_d', lb=0, ub=1e20)
            prob_nu_bar_d = model.addMVar((4, T), vtype=GRB.CONTINUOUS, name='nu_bar_d', lb=0, ub=1e20)
            prob_rou_d = model.addMVar((4, T + 1), vtype=GRB.CONTINUOUS, name='rou_d', lb=-1e20, ub=1e20)
            prob_delta_d = model.addMVar((7, T), vtype=GRB.CONTINUOUS, name='delta_d', lb=-1e20, ub=1e20)

            model.setObjective(-0.25 * sum(prob_rou_d[:, t] @ prob_rou_d[:, t] for t in range(T + 1))
                               - sum(0.25 * (prob_delta_d[:, t] @ prob_delta_d[:, t]) + prob_mu_d[:, t] @ prob_h
                                     + prob_nu_bar_d[:, t] @ prob_v_bar[:, t]
                                     - prob_nu_underline_d[:, t] @ prob_v_underline[:, t]
                                     for t in range(T)) - prob_lambda_d[:, 0] @ prob_x0, GRB.MAXIMIZE)

            for t in range(T):
                model.addConstr(prob_QT @ prob_rou_d[:, t] + prob_lambda_d[:, t] - prob_AT @ prob_lambda_d[:, t + 1]
                                + prob_FT @ prob_mu_d[:, t] == 0, name='d1')

            model.addConstr(prob_QT @ prob_rou_d[:, T] + prob_lambda_d[:, T] == 0, name='d2')

            for t in range(T):
                model.addConstr(prob_RT @ prob_delta_d[:, t] - prob_BT @ prob_lambda_d[:, t + 1]
                                + prob_GT @ prob_mu_d[:, t] + prob_VT @ prob_nu_bar_d[:, t]
                                - prob_VT @ prob_nu_underline_d[:, t] == 0, name='d3')

        return model

    def isIntegral(self, model):
        # 参数
        T = 2  # time horizon
        u_dimension = 7  # input dimension
        u_binary_dimension = 4
        # prob.optimize() 求解对偶问题，怎么知道原变量？判断是不是整数解？
        # 判断整数与否针对的是当前求解的特定的问题中的整数变量，因此因当先完成 buildProblem 函数
        # 获取 shadow price，Primal问题原变量的最优值
        d = model.getConstrs()
        dualArray = []
        horizon_is_1 = self.horizon_is_1

        # print(model.getAttr(GRB.Attr.NumConstrs))
        for i in range(model.getAttr(GRB.Attr.NumConstrs)):
            dualArray.append(d[i].getAttr(GRB.Attr.Pi))
        # print('pi:', dualArray)

        if horizon_is_1: # T = 1
            # 1. 从 price 中获取离散变量对应的值
            v = dualArray[7:11] # 因为第 7 个到第 10 个为整数变量，问题不同，参数要改
            self.bool_vars = list(map(abs, v))  # 因为是等式约束，所以有可能为正，有可能为负
        else: # T > 1，将矩阵按列合并成一个数列，存入 bool_vars
            v = np.zeros(shape=(u_binary_dimension, T))
            v_temp = dualArray[-u_dimension:]
            v[:, -1] = v_temp[-u_binary_dimension:]
            for i in range(1, T):
                v_temp = dualArray[-(i + 1) * u_dimension:-i * u_dimension]
                v[:, -i-1] = v_temp[-u_binary_dimension:]
            # 拓展时域 T 之后，bool_vars 应该是一个矩阵
            self.bool_vars = np.concatenate(np.absolute(v).T, axis= 0)

        # print('bool_vars = ',self.bool_vars)
        return all([abs(v- 1) <= 1e-3 or abs(v - 0) <= 1e-3 for v in self.bool_vars])

    def branch(self, model):
        # 循环中已经判断了该问题的最优解不是整数，所以才会进入分枝
        # 从 v_{t|\tau} 中找最小的 t，满足其中有元素不是整数，按照这个元素分成两个集合
        # 根据两个集合构造相应的节点，需要 变量，01变量，目标函数，约束，上界，下界
        # 此时传进来的 model 是优化结束的 model，要拿到它的 solution，为 children 构造初始解， 同时构造下界

        # 1. self.bool_vars 中存的是一个矩阵，列表形式实现 [[v_{0 | \tau}], [v_{1 | \tau}], ..., v_{T-1 | \tau}]
        # 找到具体要分枝的变量，列表形式实现 v_{t | \tau}， v_{t | \tau} 中有 mu 个元素，确定是哪一个，记作tbb。
        # 2. 父节点中下界 v_undeline 中 tbb 置 1，与父节点中原上界配合成为新节点
        # 3. 父节点中上界 v_bar 中 tbb 置 1， 与父节点中原下界配合成为新节点
        # 4. 计算新节点的 下界。
        bool_vars = np.array(self.bool_vars)
        horizon_is_1 = self.horizon_is_1
        # bool_vars = np.array([0.3, 0.19999999999999996, 0.09999999999999998, 0.0])
        # 将非常接近 0 和 1 的数全部置为 0
        bool_vars[bool_vars<=1e-3] = 0
        bool_vars[bool_vars>=1-1e-3] = 0
        non_integral=np.nonzero(bool_vars)

        children = []

        non_integral_d = non_integral[0][0]
        node1 = copy.deepcopy(self)
        node1.v_bar[non_integral_d] = 0
        # 下界还没求
        node1.LB = 0
        node2 = copy.deepcopy(self)
        node2.v_underline[non_integral_d] = 1
        # 下界还没求
        node2.LB = 0
        children.append(node1)
        children.append(node2)
        print('node1', node1.v_underline, node1.v_bar)
        print('node2', node2.v_underline, node2.v_bar)

        return children

# 准备动手写滚动时域的代码，
class Hybrid_MPC():
    def __init__(self, maxtime=int(), time_horizon=int(), u_binary_dimension=int()):
        self.maxtime = maxtime
        self.u_binary_dimension = u_binary_dimension
        self.time_horizon = time_horizon

    def hybridMPCSolve(self):
        # initialization
        t = 0
        time_used = []
        T = self.time_horizon
        initial_v_underline = np.zeros(shape=(self.time_horizon * self.u_binary_dimension))
        initial_v_bar = np.ones(shape=(self.time_horizon * self.u_binary_dimension))
        rootnode = BBTreeNode(v_underline=initial_v_underline, v_bar=initial_v_bar, horizon_is_1=False)
        frontier = []
        frontier.append(rootnode)
        tree = BBTree(frontier=frontier)
        u_continuous_best = []
        u_binary_best = np.zeros(shape=(self.u_binary_dimension, 1))
        while t < self.maxtime:
            # 计算时间
            time_start = time.time()

            tree.bbtreeSolve()
            # 得到现在求解的 frontier
            frontier_old = tree.frontier
            # 得到现在求解的 v_{0|\tau}
            v0 = tree.bestnode.bool_vars[0: self.u_binary_dimension]
            # 保留 tree 需要留下来的数据

            u_binary_best = np.c_[u_binary_best, v0.T]

            # 对 frontier 中的每个集合处理，构造新的 frontier
            frontier_new = []
            for node in frontier_old:
                v0_underline = node.v_underline[0: self.u_binary_dimension]
                v_underline_left = node.v_underline[self.u_binary_dimension : T * self.u_binary_dimension]
                v0_bar = node.v_bar[0: self.u_binary_dimension]
                v_bar_left = node.v_bar[self.u_binary_dimension: T * self.u_binary_dimension]

                if all(v0_underline <= v0) and all(v0 <= v0_bar): # 节点中包含了 v0
                    v_underline_new = np.append(v_underline_left, np.zeros(shape=(1, self.u_binary_dimension)))
                    v_bar_new = np.append(v_bar_left, np.ones(shape=(1, self.u_binary_dimension)))
                    # 这里创建了新的节点，如果只是改动旧节点，会不会速度更快
                    node_new = BBTreeNode(v_underline=v_underline_new, v_bar=v_bar_new, horizon_is_1=False)
                    frontier_new.append(node_new)
            # 用新的 frontier 构造 BBTree
            tree = BBTree(frontier=frontier_new)
            t += 1

            # 计算时间
            time_end = time.time()
            time_used.append(time_end-time_start)

            print('WoW! 求完第', t, '秒了！')
        return time_used


# 主代码
if __name__ == '__main__':
    controller = Hybrid_MPC(maxtime=5, time_horizon=2, u_binary_dimension=4)
    computation_time = controller.hybridMPCSolve()

    print(1)

# frontier = []
# node = BBTreeNode(v_underline=np.array([0,0,0,0,0,0,0,0]), v_bar=np.array([1,1,1,1,1,1,1,1]), horizon_is_1=False)
# frontier.append(node)
# tree = BBTree(frontier=frontier)
# tree.bbtreeSolve()
# print(tree.bestnode.LB)
# print(tree.bestnode.UB)
# print(tree.bestnode.bool_vars)

# 没有正式做拉格朗日的下界，目前求解这个需要搜索45个节点，256个需要搜索45个，如果是把下界加进来呢？
# 会更快剪枝吗？

# prob = node.buildProblem()
# prob.optimize()
# node.isIntegral(prob)
# print('obj=', prob.objVal)
# print('is_Integral=', node.isIntegral(prob))
#
# for v in prob.getVars():
#     print(v.varName, ':', v.x)














































