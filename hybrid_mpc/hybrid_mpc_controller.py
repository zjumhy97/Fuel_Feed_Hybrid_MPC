import gurobipy as gp
from gurobipy import GRB
import numpy as np

__all__ = ['HybridMPCController']

class HybridMPCController:
    def __init__(self, discrete_dim, continuous_dim, control_horizon, prediction_horizon) -> None:
        self.discrete_dim = discrete_dim
        self.continuous_dim = continuous_dim
        self.var_dim = discrete_dim + continuous_dim
        self.control_horizon = control_horizon
        self.prediction_horizon = prediction_horizon
        if self.control_horizon > self.prediction_horizon:
            raise ValueError("The control horizon must be leq to the prediction horizon!")
        
        # initialize the timer
        self.timer = [0] * self.var_dim
        

    def setup(self):
        model = gp.Model("HMPC")

        # Create variables
        # x = m.addVar(vtype=GRB.BINARY, name="x") 
        # y = m.addVar(vtype=GRB.CONTINUOUS, name="y")

        # model.setObjective(, GRB.MINIMIZE)

        # model.addConstr(, "c0")
        # model.addConstr(, "c1")

        # model.setObjective(, GRB.MINIMIZE)

        model.optimize()

        for v in m.getVars():
            print('%s %g' % (v.varName, v.x))

        print('Obj: %g' % m.objVal)


    def _solve(self, state: np.array) -> np.array:
        """
        Solve the optimization problem in each circle.
        
        Input:
        state: the state of the system. 
        
        Output:
        action: the action will be applied on the system.
        """
        pass
        return action


    def solve(self):
        """
        The algorithm of the Hybrid MPC.
        """
        pass

if __name__ == '__main__':
    controller_param = {
        'discrete_dim': 6,
        'continuous_dim': 6,
        'control_horizon': 1,
        'prediction_horizon': 1
    }
    controller = HybridMPCController(**controller_param)
    print(1)
