import gurobipy.gurobipy as gp
from gurobipy import *
from gurobipy.gurobipy import GRB
import numpy as np
from params import configs

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

        self.output_signal = []
        
    def setup(self, A: np.ndarray, B: np.ndarray, x_tau: np.ndarray, x_ub:np.ndarray, u_ub:np.ndarray, h:np.ndarray, ideal_cg:np.ndarray, M_hat:float, cg_of_tanks:np.ndarray, timer:np.ndarray, fuel_feed_system) -> np.ndarray:
        model = gp.Model("HMPC")

        # system state
        x = model.addMVar((self.system_state_dim, self.prediction_horizon+1), vtype=GRB.CONTINUOUS, name="x", lb=0)

        # controller input
        u_binary = model.addMVar((self.discrete_dim, self.prediction_horizon), vtype=GRB.BINARY, name="u_binary")
        u_continuous = model.addMVar((self.continuous_dim, self.prediction_horizon), vtype=GRB.CONTINUOUS, name="u_continuous", lb=0)

        # dynamics constraint
        model.addConstr(x[:, 0] == x_tau, "c:x_0")
        for t in range(self.prediction_horizon):
            model.addConstr(x[:, t+1] == A @ x[:, t] + B @ u_continuous[:, t], name="c:x_t")

        # state range constraint
        model.addConstrs(x[:, t] <= x_ub for t in range(self.prediction_horizon+1))

        # continuous input upper bound
        model.addConstrs(u_continuous[:, t] <= u_ub for t in range(self.prediction_horizon))
        K = 1000
        model.addConstrs(u_continuous[:, t] <= K * u_binary[:, t] for t in range(self.prediction_horizon))
        
        # test:
        for i in range(self.discrete_dim):
            if fuel_feed_system.tanks[i+1].mass == 0:
                model.addConstrs(u_binary[i, t] == 0 for t in range(self.prediction_horizon))
                model.addConstrs(u_continuous[i, t] == 0 for t in range(self.prediction_horizon))

        # fuel feed structural bound
        model.addConstrs(gp.quicksum(u_binary[i, t] for i in range(2, 6)) <= 2 for t in range(self.prediction_horizon))
        model.addConstrs(gp.quicksum(u_binary[:, t]) <= 3 for t in range(self.prediction_horizon))

        # engine fuel feed constraint
        # attention: not all the tanks!
        model.addConstrs(gp.quicksum(u_continuous[j, t] for j in [1, 2, 3, 4]) >= h[t] for t in range(self.prediction_horizon))
        model.addConstrs(gp.quicksum(u_continuous[j, t] for j in [1, 2, 3, 4]) <= h[t] * (1.05) for t in range(self.prediction_horizon))

        # constraints for timer
        for i in range(self.discrete_dim):
            if 0 < timer[i] < configs.minimum_duration:
                model.addConstrs(u_binary[i, t] == 1 for t in range(min(self.prediction_horizon, int(configs.minimum_duration - timer[i]))))

        # todo: 
        objective = 0
        for t in range(self.prediction_horizon+1):
            # predict_cg = 1 / M_hat * (x[:,t].reshape(1,-1) @ cg_of_tanks)[0]
            predict_cg =  (x[:,t].reshape(1,-1) @ cg_of_tanks)[0]
            error = predict_cg - ideal_cg[t] * M_hat
            objective += gp.quicksum(error[i]*error[i] for i in range(3))

        model.setObjective(objective, GRB.MINIMIZE)
        model.setParam('OutputFlag', configs.output_flag)
        model.optimize()

        # print(x.x)
        # print(u_binary.x)
        # print(u_continuous.x)
        # for v in m.getVars():
        #     print('%s %g' % (v.varName, v.x))

        # print('Obj: %g' % m.objVal)
        u_binary_output = u_binary.x[:, 0]
        for i in range(configs.discrete_dim):
            u_binary_output[i] = int(u_binary_output[i])
        u_continuous_output = u_continuous.x[:, 0]
        output_signal = np.append(u_binary_output, u_continuous_output)
        self.output_signal.append(list(output_signal))
        
        # print('The output:', output_signal)
        return output_signal
      

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
