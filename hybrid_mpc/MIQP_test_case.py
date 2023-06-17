import numpy as np
import gurobipy as gp
from gurobipy import *

try:
    # 参数
    Q = np.identity(4) # Q^TQ=Q
    x0 = np.array([0.3, 0.4, 0.5, 0.6])
    A = np.identity(4)
    B = np.array([[0.1, 0.2, 0.3, 1, 0, 0, 0],
                  [0.1, 0.2, 0.3, 0, 1, 0, 0],
                  [0.1, 0.2, 0.3, 0, 0, 1, 0],
                  [0.1, 0.2, 0.3, 0, 0, 0, 1]])
    BT = np.transpose(B)
    F = np.identity(4)
    G = np.array([[0, 0, 0, 1, 0, 0, 0],
                  [0, 0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 1]])
    GT = np.transpose(G)
    V = np.array([[0, 0, 0, 1, 0, 0, 0],
                  [0, 0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 1]])
    VT = np.transpose(V)
    v_bar = np.array([1, 1, 1, 1])
    v_underline = np.array([0, 0, 0, 0])
    h = np.array([1, 1, 1.5, 1])

#     # 模型
#     model = gp.Model('miqp')
#
#     # 变量
#     x = model.addMVar((4, 2), vtype=GRB.CONTINUOUS, name='x',lb=-1e20,ub=1e20)
#     u = model.addMVar(7, vtype=GRB.CONTINUOUS, name='u',lb=-1e20,ub=1e20)
#
#     # 目标函数
#     # model.setObjective(x[:,0] @ Q @ x[:,0] + x[:,1] @ Q @ x[:,1], GRB.MINIMIZE)
#     model.setObjective(sum(x[:, i] @ Q @ x[:, i] for i in range(2)), GRB.MINIMIZE)
#     # 约束
#     model.addConstr(x[:, 0] == x0, name='c1')
#     model.addConstr(x[:, 1] == A @ x[:, 0] + B @ u, name='c2')
#     model.addConstr(F @ x[:, 0] + G @ u <= h, name='c3')
#     model.addConstr(V @ u >= v_underline, name='c4')
#     model.addConstr(V @ u <= v_bar, name='c5')
#     # 求解
#     # model.setParam('outPutFlag', 0)  # 不输出求解日志
#     model.optimize()
#
#     # 输出
#     print('obj=', model.objVal)
#     for v in model.getVars():
#         print(v.varName, ':', v.x)
#
# except GurobiError as e:
#     print('Error code ' + str(e.errno) + ':' + str(e))
#
# except AttributeError:
#     print('Encountered an attribute error')


    # the optimum should be 1.04750001
    # Dual part
    # 模型
    model = gp.Model('miqp')

    # 变量
    # lambda_d = model.addVars((4, 2), vtype=GRB.CONTINUOUS, name='lambda_d')
    # mu_d = model.addVars(4, vtype=GRB.CONTINUOUS, name='mu_d', lb=0)
    # nu_d = model.addVars((4, 2), vtype=GRB.CONTINUOUS, name='nu_d', lb=0)

    lambda_d = model.addMVar((4, 2), vtype=GRB.CONTINUOUS, name='lambda_d',lb=-1e20,ub=1e20)
    mu_d = model.addMVar(4, vtype=GRB.CONTINUOUS, name='mu_d', lb=0, ub=1000)
    nu_d = model.addMVar((4, 2), vtype=GRB.CONTINUOUS, name='nu_d', lb=0,ub=1e20)
    z_d = model.addMVar(4, vtype=GRB.CONTINUOUS, name='lambda_d',lb=-1e20,ub=1e20)

    # 目标函数
    # model.setObjective(- 0.25 * ((lambda_d[:, 0] - A @ lambda_d[:, 1] + F @ mu_d) @
    #                    (lambda_d[:, 0] - A @ lambda_d[:, 1] + F @ mu_d))
    #                    - 0.25 * (lambda_d[:, 1] @ lambda_d[:, 1]) - lambda_d[:, 0] @ x0
    #                    - mu_d @ h + nu_d[:, 0] @ v_underline - nu_d[:, 1] @ v_bar, GRB.MAXIMIZE)

    # model.setObjective((lambda_d[:, 0] - A @ lambda_d[:, 1]) @ (lambda_d[:, 0] - A @ lambda_d[:, 1])
    #                    , GRB.MINIMIZE)

    # model.setObjective(-0.25 * (lambda_d[:, 0] @ lambda_d[:, 0] + ), GRB.MINIMIZE) # 这个没问题

    model.setObjective(-0.25 * (z_d @ z_d)
                       - 0.25 * (lambda_d[:, 1] @ lambda_d[:, 1]) - lambda_d[:, 0] @ x0
                       - mu_d @ h + nu_d[:, 0] @ v_underline - nu_d[:, 1] @ v_bar
                       , GRB.MAXIMIZE)

    # 约束
    model.addConstr(z_d == lambda_d[:, 0] - A @ lambda_d[:, 1] + F @ mu_d, name='d1')
    model.addConstr(-BT @ lambda_d[:, 1] + GT @ mu_d + VT @ nu_d[:, 1] - VT @ nu_d[:, 0] == 0, name='d2')

    # 求解
    # model.setParam('outPutFlag', 0)  # 不输出求解日志
    model.optimize()

    # 输出
    print('obj=', model.objVal)
    for v in model.getVars():
        print(v.varName, ':', v.x)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ':' + str(e))

except AttributeError:
    print('Encountered an attribute error')

