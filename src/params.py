import argparse

parser = argparse.ArgumentParser(description='Arguments for Fuel_Feed_Hybrid_MPC')

parser.add_argument('--oil_density', type=float, default=850.0, help='The density of oil.')
parser.add_argument('--aircraft_net_mass', type=float, default=3000.0, help='The net mass of the aircraft.')

# controller
parser.add_argument('--discrete_dim', type=int, default=6, help='The discrete dimension of the controller.')
parser.add_argument('--continuous_dim', type=int, default=6, help='The continuous dimension of the controller.')
parser.add_argument('--control_horizon', type=int, default=1, help='The control horizon of the controller.')
parser.add_argument('--prediction_horizon', type=int, default=5, help='The prediction horizon of the controller.')
parser.add_argument('--system_state_dim', type=int, default=6, help='The system state dimension of the controller.')
parser.add_argument('--minimum_duration', type=int, default=60, help='The minimum duration of the fuel feed.')
parser.add_argument('--time_length', type=int, default=7200, help='The time length of the controller working.')

# optimization problem
parser.add_argument('--output_flag', type=int, default=0, help='Whether print the optimization log.')

# printer
parser.add_argument('--print_cycle', type=int, default=1, help='The print circle of the printer.')

configs = parser.parse_args()