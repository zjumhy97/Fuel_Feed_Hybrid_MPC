# bb_tree_gurobi.py
# author: Haoyu Miao
# date: 2022/02/24

# import gurobipy

import cvxpy as cvx
import gurobipy as gp
import copy
from heapq import *
import numpy as np
import itertools
counter = itertools.count() # create a iterator from 0

class BBTreeNode():
    def __init__(self, vars=set(), constraints=[], objective=0, bool_vars=set()):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.bool_vars = bool_vars
        self.children = [] # Q: what is the use of children?
    def buildProblem(self):
        prob = cvx.Problem(cvx.Minimize(self.objective), self.constraints)
        #i put Minimize, just so you know that I'm assuming it
        return prob
    def is_integral(self):
        # Q: why need to determine whether the variable is integer?
        return all([abs(v.value - 1) <= 1e-3 or abs(v.value - 0) <= 1e-3 for v in self.bool_vars])
    def branch(self):
        children = []
        for b in [0,1]:
            n1 = copy.deepcopy(self) #yeesh. Not good performance wise, but is simple implementation-wise
            v = n1.heuristic() #dangerous what if they don't do the same one?
            # I need to do it here though because I need access to copied v.
            n1.constraints.append( v == b ) # add in the new binary constraint
            n1.children = []
            n1.bool_vars.remove(v) #remove binary constraint from bool var set
            n1.vars.add(v) #and add it into var set for later inspection of answer
            #self.children.append(n1)   # eventually I might want to keep around the entire search tree.
            # I messed this up though
            children.append(n1)
        return children
    def heuristic(self):
        # a basic heuristic of taking the ones it seems pretty sure about
        return min([(min(1 - v.value, v.value) , i, v) for i, v in enumerate(self.bool_vars)])[2]
    def bbsolve(self):
        root = self
        # solve the B&B tree root node problem (relaxation)
        # 作者这里是为了求解第一个节点，即根节点
        res = root.buildProblem().solve()
        # 根节点压入heap
        heap = [(res, next(counter), root)]
        bestres = 1e20 # a big arbitrary initial best objective value ??? self.treeUB
        bestnode = root # initialize bestnode to the root ???
        print(heap)
        nodecount = 0
        while len(heap) > 0:
            nodecount += 1 # for statistics
            print("Heap Size: ", len(heap))
            _, _, node = heappop(heap) # Q: what is the use of heappop?
            prob = node.buildProblem()
            res = prob.solve()
            print("Result: ", res)
            if prob.status not in ["infeasible", "unbounded"]:
                # Q: if the problem state is "infeasible", what operation should be done?
                # Q: if the problem state is "unbounded", what operation should be done?
                if res > bestres - 1e-3: #even the relaxed problem sucks. forget about this branch then
                    print("Relaxed Problem Stinks. Killing this branch.")
                    pass
                elif node.is_integral(): #if a valid solution then this is the new best
                        print("New Best Integral solution.")
                        bestres = res
                        bestnode = node
                else:
                    #otherwise, we're unsure if this branch holds promise.
                    # Maybe it can't actually achieve this lower bound. So branch into it
                    new_nodes = node.branch()
                    for new_node in new_nodes:
                        # using counter to avoid possible comparisons between nodes. It tie breaks
                        heappush(heap, (res, next(counter), new_node))
        print("Nodes searched: ", nodecount)
        return bestres, bestnode




# a toy example with the dimension of the integer variable is 3
N = 3
# prices = -np.random.rand(N)
prices = [-1, -1, 1]
# sizes = np.random.rand(N)
# print(prices)
x = cvx.Variable(N)
constraints = []
constraints += [x <= 1, 0 <= x] #The relaxation of the binary variable constraint
# constraints += [sizes*x <= 5] # total size of knapsack is 5
objective = prices * x
# objective = cvx.norm(x)
bool_vars = {x[i] for i in range(N)}
root = BBTreeNode(constraints = constraints, objective= objective, bool_vars = bool_vars)
res, sol = root.bbsolve()
print(sorted(list([(v.name(), v.value) for v in sol.bool_vars] + [(v.name(), v.value) for v in sol.vars] ) ))
























