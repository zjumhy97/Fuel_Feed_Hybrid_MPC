import gurobipy as gp
from gurobipy import GRB
import numpy as np

__all__ = ['HybridMPCController']

class HybridMPCController:
    def __init__(self, discrete_dim, continuous_dim, control_horizon, prediction_horizon, system_state_dim, ) -> None:
        self.discrete_dim = discrete_dim
        self.continuous_dim = continuous_dim
        self.var_dim = discrete_dim + continuous_dim
        self.control_horizon = control_horizon
        self.prediction_horizon = prediction_horizon
        self.system_state_dim = system_state_dim
        if self.control_horizon > self.prediction_horizon:
            raise ValueError("The control horizon must be leq to the prediction horizon!")
        
        # initialize the timer
        self.timer = [0] * self.var_dim
        

    def setup(self, A, B, x_tau: np.array):
        model = gp.Model("HMPC")

        # system state
        x = model.addMVar((self.system_state_dim, self.prediction_horizon+1), vtype=GRB.CONTINUOUS, name="x", lb=0)
        x_ub = np.array()

        # controller input
        u_binary = model.addMVar((self.discrete_dim, self.prediction_horizon), vtype=GRB.BINARY, name="u_binary")
        u_continuous = model.adMVar((self.continuous_dim, self.prediction_horizon), vtype=GRB.CONTINUOUS, name="u_continuous", lb=0)
        u_ub = np.array()

        # dynamics constraint
        model.addConstr(x[:, 0] == x_tau, "x_0")
        for t in range(self.prediction_horizon):
            model.addConstr(x[:, t+1] == A @ x[:, t] + B @ u_continuous[:, t], "x_t")

        # state range constraint
        model.addConstrs(x[:, t] <= x_ub for t in range(self.prediction_horizon+1))

        # continuous input upped bound
        model.addConstrs(u_continuous[:, t] <= u_ub for t in range(self.prediction_horizon))



        model.setObjective(0, GRB.MINIMIZE)

        model.optimize()

        # for v in m.getVars():
        #     print('%s %g' % (v.varName, v.x))

        # print('Obj: %g' % m.objVal)


    def _solve(self, state: np.array) -> np.array:
        """
        Solve the optimization problem in each circle.
        
        Input:
        state: the state of the system. 
        
        Output:
        action: the action will be applied on the system.
        """
        pass
        # return action


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
        'prediction_horizon': 1,
        'system_state_dim': 6
    }
    controller = HybridMPCController(**controller_param)
    print(1)
