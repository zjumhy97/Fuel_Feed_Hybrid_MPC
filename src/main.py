import numpy as np
import json
import pandas as pd
import time
from openpyxl import Workbook

from hybrid_mpc import HybridMPCController
from fuel_feed_system import FuelFeedSystem
from params import configs

def solve_MPC(time_length:int, controller:HybridMPCController, fuel_feed_system:FuelFeedSystem):
    results = []
    timer = np.zeros(configs.discrete_dim)
    A = np.eye(configs.system_state_dim)
    B = fuel_feed_system.topology_matrix
    x_ub = fuel_feed_system.fuel_mass_upper_bound
    u_ub = fuel_feed_system.fuel_feed_velocity_upper_bound

    for tau in range(time_length):
        time_0 = time.time()
        x_tau = fuel_feed_system.fuel_mass()
        h = fuel_consume_velocity[tau:tau+configs.prediction_horizon]
        M_hat = fuel_feed_system.total_fuel_mass() + configs.aircraft_net_mass
        # most time cost
        time_1 = time.time()
        cg_of_tanks = fuel_feed_system.center_of_gravity_of_tanks()
        time_2 = time.time()
        action = controller.setup(A=A, B=B, x_tau=x_tau, x_ub=x_ub, u_ub=u_ub, h=h, ideal_cg=ideal_cg[tau:tau+controller.prediction_horizon+1], M_hat=M_hat, cg_of_tanks=cg_of_tanks, timer=timer, fuel_feed_system=fuel_feed_system)
        fuel_feed_system.step(fuel_feed=action[controller.discrete_dim:])
        time_3 = time.time()
        timer = update_timer(timer=timer, action=action)
        # record the data
        center_of_gravity = fuel_feed_system.center_of_gravity()
        cg_difference = np.linalg.norm(np.array(fuel_feed_system.center_of_gravity() - ideal_cg[tau], dtype=np.float32))
        # results[tau] = {
        #     "action": list(action),
        #     "center_of_gravity": list(fuel_feed_system.center_of_gravity()),
        #     "cg_difference": cg_difference
        # }
        results.append([tau, action.tolist(), center_of_gravity.tolist(),float(cg_difference)])
        if tau % configs.print_cycle == 0:
            print(">> ", tau, ":", time_3 - time_0)
            print("-- timer:", timer)
            print("-- ", action)
            print("-- ", cg_difference)
            if cg_difference is None:
                print(1)
            print("\n")
        
    return results


def update_timer(timer:np.ndarray, action:np.ndarray) -> np.ndarray:
    assert len(timer) == configs.discrete_dim
    for i in range(configs.discrete_dim):
        timer[i] += int(action[i])
        if timer[i] >= configs.minimum_duration:
            timer[i] = 0
    return timer


if __name__ =='__main__':
    controller_param = {
        'discrete_dim': configs.discrete_dim,
        'continuous_dim': configs.continuous_dim,
        'control_horizon': configs.control_horizon,
        'prediction_horizon': configs.prediction_horizon,
        'system_state_dim': configs.system_state_dim
    }
    controller = HybridMPCController(**controller_param)
    
    fuel_tank_file_name = './data/fuel_tank_params.json'
    with open(fuel_tank_file_name, 'r') as file:
        json_data = json.load(file)
    fuel_feed_system = FuelFeedSystem(json_data=json_data)

    # load flight data
    flight_file_name = './data/附件3-问题2数据.xlsx'
    df1 = pd.read_excel(flight_file_name, sheet_name="发动机耗油速度")
    fuel_consume_velocity = np.array(df1['耗油速度(kg/s)'])
    df2 = pd.read_excel(flight_file_name, sheet_name='飞行器理想质心数据')
    ideal_cg = np.array(df2.loc[:,"X坐标（米）":"z坐标（米）"])

    # 7200 - 4
    # time_length = len(df1) - configs.prediction_horizon - 1
    time_length = configs.time_length

    results = solve_MPC(time_length=time_length, controller=controller, fuel_feed_system=fuel_feed_system)
    # controller.setup(A=A, B=B, x_tau=x_tau, x_ub=x_ub, u_ub=u_ub, h=h, ideal_cg=ideal_cg[0:11], M_hat=M_hat, cg_of_tanks=cg_of_tanks)
    
    f=open("./data/results.txt","w") 
    for line in results:
        f.write(str(line)+'\n')
    f.close()

    print(1)