from hybrid_mpc.hybrid_mpc_controller import HybridMPCController
from fuel_feed_system import FuelFeedSystem

fuel_feed_system_param = {
    'aircraft_weight': 10000.0
}


controller_param = {
    'discrete_dim': 6,
    'continuous_dim': 6,
    'control_horizon': 1,
    'prediction_horizon': 1
}


if __name__ =='__main__':
    controller = HybridMPCController(**controller_param)
    fuel_feed_system = FuelFeedSystem(**fuel_feed_system_param)