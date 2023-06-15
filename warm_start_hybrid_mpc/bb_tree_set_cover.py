# bb_tree_set_cover.py
# author: Haoyu Miao
# date: 2022/03/02

import gurobipy as gp
from gurobipy import GRB
from operator import attrgetter
import copy
from heapq import *
import numpy as np
import itertools
counter = itertools.count() # create a iterator from 0

class BBTree():
    def __init__(self, root_nodes=set(), tolerance=1e-3, treeUB=float('inf'), treeLB=-float('inf')):
        # 按照 paper 中的说法，frontier，其中每个集合的下界，以及 treeUB，这三个是必须的
        # 【重要！！！】 frontier 中每个 node 的下界应该是 parent node 对应对偶问题的界
        self.frontier = root_nodes
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
            heappush(heap, node)
        # print(heap)
        # Q: 根节点需要求解吗？
        # Q: 需要更新B&B树上界和最新节点吗？
        nodecount = 0
        while len(heap) > 0: # 这个循环是边界 frontier 迭代的循环，第0次，frontier就是根节点
            print('-----------  ', nodecount, '  -----------')
            # 用 len(heap) 作为循环条件的缺点——每次求解bb树都要求解完，而不是中途求到feasible的解也可以退出
            nodecount += 1
            print("Heap Size: ", len(heap))
            # solve the subproblem chosen by selectSubproblem function, get relaxation value \theta(V^{i})
            # 满足条件的话，从边界中选一个节点求解, 2:optimal 3:infeasible 5:unbounded
            node = self.selectSubproblem(heap) # 从 frontier 中还是 heap 中？
            # 选择一个节点后，应该把节点从 heap 中弹出
            heap.remove(node) # heap 中移除，但是 frontier 中不移除

            # buildProblem 这里需要参数的哦~
            prob = node.buildProblem()
            prob.optimize()
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
                else: # branching
                    # 现在开始处理分支，争取2022/03/05 23:30前处理完
                    newnodes = node.branch(model=prob)
                    self.frontier.remove(node)
                    for newnode in newnodes:
                        self.frontier.append(newnode)
                        heap.append(newnode)
                        # heappush(self.frontier, newnode)
                        # heappush(heap, newnode)
        print("Nodes searched: ", nodecount)
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
        prob_Q = np.identity(4)  # Q^TQ=Q
        # 初始状态需要作为参数吗？因为涉及到每个时间步时，具体的初始状态都是不一样的。
        prob_x0 = np.array([0.3, 0.4, 0.5, 0.6])
        # prob_x0 = x0
        prob_A = np.identity(4)
        prob_B = np.array([[0.1, 0.2, 0.3, 1, 0, 0, 0],
                           [0.1, 0.2, 0.3, 0, 1, 0, 0],
                           [0.1, 0.2, 0.3, 0, 0, 1, 0],
                           [0.1, 0.2, 0.3, 0, 0, 0, 1]])
        prob_BT = np.transpose(prob_B)
        prob_F = np.identity(4)
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
        prob_h = np.array([1, 1, 1.5, 1])
        # 二元变量线性松弛的上下界
        # 在 B&B 框架中，最重要的应该是下面这两个变量，这两个变量是二元变量的上下界，这两个应该是传进来的
        # prob_v_bar = np.array([1, 1, 1, 1])
        # prob_v_underline = np.array([0, 0, 0, 0])
        prob_v_bar = self.v_bar
        prob_v_underline = self.v_underline

        model = gp.Model('miqp_dual')
        model.setParam('Method', 1)
        model.setParam('outPutFlag', 0)

        prob_lambda_d = model.addMVar((4, 2), vtype=GRB.CONTINUOUS, name='lambda_d', lb=-1e20, ub=1e20)
        prob_mu_d = model.addMVar(4, vtype=GRB.CONTINUOUS, name='mu_d', lb=0, ub=1000)
        prob_nu_d = model.addMVar((4, 2), vtype=GRB.CONTINUOUS, name='nu_d', lb=0, ub=1e20)
        prob_z_d = model.addMVar(4, vtype=GRB.CONTINUOUS, name='z_d', lb=-1e20, ub=1e20)

        model.setObjective(-0.25 * (prob_z_d @ prob_z_d)
                           - 0.25 * (prob_lambda_d[:, 1] @ prob_lambda_d[:, 1]) - prob_lambda_d[:, 0] @ prob_x0
                           - prob_mu_d @ prob_h + prob_nu_d[:, 0] @ prob_v_underline - prob_nu_d[:, 1] @ prob_v_bar
                           , GRB.MAXIMIZE)
        # 变量替换约束（方便gurobipy编程）
        model.addConstr(prob_z_d == prob_lambda_d[:, 0] - prob_A @ prob_lambda_d[:, 1] + prob_F @ prob_mu_d, name='d1')
        # 最优性约束
        model.addConstr(
            -prob_BT @ prob_lambda_d[:, 1] + prob_GT @ prob_mu_d + prob_VT @ prob_nu_d[:, 1] - prob_VT @ prob_nu_d[:,
                                                                                                    0] == 0,name='d2')
        # --------------------------------------------------------------------------------------------------------------
        # 下面这个是问题的 Primal 形式，但是这篇文章中用的是对偶问题
        # prob_v_bar = v_bar
        # prob_v_underline = v_underline
        # model = gp.Model('miqp')
        #
        # prob_x = model.addMVar((4, 2), vtype=GRB.CONTINUOUS, name='x', lb=-1e20, ub=1e20)
        # prob_u = model.addMVar(7, vtype=GRB.CONTINUOUS, name='u', lb=-1e20, ub=1e20)
        #
        # model.setObjective(sum(prob_x[:, i] @ prob_Q @ prob_x[:, i] for i in range(2)), GRB.MINIMIZE)
        #
        # model.addConstr(prob_x[:, 0] == prob_x0, name='c1')
        # model.addConstr(prob_x[:, 1] == prob_A @ prob_x[:, 0] + prob_B @ prob_u, name='c2')
        # model.addConstr(prob_F @ prob_x[:, 0] + prob_G @ prob_u <= prob_h, name='c3')
        # model.addConstr(prob_V @ prob_u >= prob_v_underline, name='c4')
        # model.addConstr(prob_V @ prob_u <= prob_v_bar, name='c5')
        # ---------------------------------------------------------------------------------------------------------------
        return model

    def isIntegral(self, model):
        # prob.optimize() 求解对偶问题，怎么知道原变量？判断是不是整数解？
        # 判断整数与否针对的是当前求解的特定的问题中的整数变量，因此因当先完成 buildProblem 函数
        # 获取 shadow price，Primal问题原变量的最优值
        d = model.getConstrs()
        dualArray = []
        print(model.getAttr(GRB.Attr.NumConstrs))
        for i in range(model.getAttr(GRB.Attr.NumConstrs)):
            dualArray.append(d[i].getAttr(GRB.Attr.Pi))
        print('pi:', dualArray)
        # 1. 从 price 中获取离散变量对应的值
        v = dualArray[7:11] # 因为第 7 个到第 10 个为整数变量
        # 拓展时域 T 之后，bool_vars 应该是一个矩阵
        self.bool_vars = list(map(abs,v)) # 因为是等式约束，所以有可能为正，有可能为负
        print('bool_vars = ',self.bool_vars)
        # 2. 判断这些值是否为整数，如果为整数，is_Integral 置为 True
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
        bool_vars[bool_vars<=1e-3] = 0
        bool_vars[bool_vars>=1-1e-3] = 0
        non_integral=np.nonzero(bool_vars)

        children = []
        if horizon_is_1: # T=1
            non_integral_d = non_integral[0][0]
            node1 = copy.deepcopy(self)
            node1.v_bar[non_integral_d] = 0
            # 下界还没求
            node1.LB = 0
            node2 = copy.deepcopy(self)
            node2.v_underline[non_integral_d] = 1
            node2.LB = 0
            children.append(node1)
            children.append(node2)
            print('node1', node1.v_underline, node1.v_bar)
            print('node2', node2.v_underline, node2.v_bar)

        else: #T>1
            non_integral_t, non_integral_d = non_integral[0][0], non_integral[1][0]

        return children





# 主代码
frontier = []
node = BBTreeNode(v_underline=np.array([0,0,0,0]), v_bar=np.array([1,1,1,1]))
frontier.append(node)
tree = BBTree(root_nodes=frontier)
tree.bbtreeSolve()
# prob = node.buildProblem(x0=np.array([0.3, 0.4, 0.5, 0.6]), )
# prob.optimize()
# print('obj=', prob.objVal)
# print('is_Integral=', node.isIntegral(prob))

# for v in prob.getVars():
#     print(v.varName, ':', v.x)




[0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0]
[0.0, 2.9999999999999996, 0.0, -0.3, -0.19999999999999996, -0.09999999999999998, -0.0]










































