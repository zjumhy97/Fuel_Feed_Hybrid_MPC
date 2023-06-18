import argparse

parser = argparse.ArgumentParser(description='Arguments for Fuel_Feed_Hybrid_MPC')

parser.add_argument('--oil_density', type=float, default=850.0, help='The density of oil.')
parser.add_argument('--aircraft_net_mass', type=float, default=3000.0, help='The net mass of the aircraft.')
# parser.add_argument('--', type=float, default=, help='')

configs = parser.parse_args()

