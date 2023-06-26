import numpy as np
import sympy as sp
from functools import cached_property
import json

from params import configs

__all__ = ['FuelFeedSystem']

class FuelFeedSystem:
    def __init__(self, file_name:str,) -> None:
        self.tanks_directly_to_engine = []
        with open(file_name, 'r') as file:
            json_data = json.load(file)
        self.num_of_tanks = len(json_data)
        self.tanks = ['engine']
        for key in json_data:
            new_tank = FuelTank(id_number=int(key),**json_data[key])
            self.tanks.append(new_tank)

    @cached_property
    def topology_matrix(self) -> np.ndarray:
        """
        Obtain the topology matrix of the fuel feed system.
        """
        topology_matrix = np.zeros([self.num_of_tanks, self.num_of_tanks])
        for i in range(1, self.num_of_tanks + 1):
            if len(self.tanks[i].in_id) > 0:
                j = self.tanks[i].in_id[0]
                topology_matrix[i-1][j-1] = 1
            if len(self.tanks[i].out_id) > 0:
                j = self.tanks[i].out_id[0]
                topology_matrix[i-1][i-1] = -1
                if j > 0:
                    topology_matrix[j-1][i-1] = 1
                elif j == 0:
                    self.tanks_directly_to_engine.append(i)
        return topology_matrix

    @cached_property
    def output_fuel_mass_upper_bound(self) -> np.ndarray:
        ub = np.zeros((self.num_of_tanks,1))
        for i in range(self.num_of_tanks):
            ub[i] = self.tanks[i+1].velocity_UB
        return ub

    def center_of_gravity(self) -> np.ndarray:
        """
        Calculate the CG of the total aircraft.
        """
        total_mass = configs.aircraft_net_mass
        mr = 0
        for i in range(1, self.num_of_tanks+1):
            total_mass += self.tanks[i].mass
            cg_of_tank = self.tanks[i].center_of_gravity()
            mr += self.tanks[i].mass * cg_of_tank
        return mr/total_mass

    def step(self, fuel_feed: np.ndarray) -> None:
        """
        Update the fuel mass of the tanks.
        
        Input: 
        fuel_feed: The fuel feed vector, unit: kg/s
        """
        assert len(fuel_feed) == self.num_of_tanks
        variation = np.dot(self.topology_matrix, fuel_feed)
        for i in range(1, self.num_of_tanks+1):
            self.tanks[i].mass = self.tanks[i].mass + variation[i-1]
            if self.tanks[i].mass < 0:
                raise ValueError("The fuel feeded by tank", str(i), "exceeds its capacity.")



class FuelTank:
    def __init__(self, id_number:int, shape_size:list, origin:list, velocity_UB:float, initial_oil_volume:float, in_id:list, out_id:list) -> None:
        """
        Input:
        shape_size: three-demension vector, the size of the tank
        origin: the original point of the tank's coordinate frame w.r.t flight coordinate frame
        """
        self.id = id_number
        self.shape_size = np.array(shape_size)
        self.origin = np.array(origin)
        self.velocity_UB = velocity_UB
        self.oil_volume = initial_oil_volume
        self.mass = configs.oil_density * self.oil_volume
        self.theta = 0.0
        self.in_id = in_id
        self.out_id = out_id
        self.a = shape_size[0]
        self.b = shape_size[1]
        self.c = shape_size[2]
        self.V = self.a * self.b * self.c
        self.M = configs.oil_density * self.V
    
    def center_of_gravity(self) -> np.ndarray:
        """
        Compute the center of gravity (CG) of a tank
        
        Input:
        mass: the mass of oil in the tank
        theta: the pitch angle of the flight

        Output:
        cg: the center of gravity (CG) of a tank
        """
        sign = 0 if self.theta < 0 else 1

        # case 1
        xA, zB = sp.symbols('xA zB')
        eqns = [zB - 0.5 * self.c + np.tan(self.theta) * (0.5 * self.a * sign - xA)]
        eqns.append((self.M - self.mass) / self.M - 0.5 * (0.5 * self.a - sign * xA) * (0.5 * self.c - zB) / (self.a * self.c))
        S = sp.solve(eqns, [xA, zB], dict=True, real=True)
        if len(S) == 1:
            # print('case 1')
            xA_val = S[0][xA]
            zB_val = S[0][zB]
            if xA_val >= -0.5 * self.a and xA_val <= 0.5 * self.a and zB_val >= -0.5 * self.c and zB_val <= 0.5 * self.c:
                cg = -(self.M - self.mass) / self.mass * np.array([(xA_val + self.a * sign) / 3, 0, (zB_val + self.c) / 3])
                cg[0] = cg[0] + self.origin[0]
                cg[1] = cg[1] + self.origin[1]
                cg[2] = cg[2] + self.origin[2]
                return cg

        # case 2
        xB = sp.symbols('xB')
        eqns = [np.tan(self.theta) * (xA - xB) + self.c]
        eqns.append(self.mass / self.M - 0.5 - 0.5/self.a*(xA+xB)*sign)
        S = sp.solve(eqns, [xA, xB], dict=True, real=True)
        if len(S) == 1:
            # print('case 2')
            xA_val = S[0][xA]
            xB_val = S[0][xB]
            if xA_val >= -0.5 * self.a and xA_val <= 0.5 * self.a and xB_val >= -0.5 * self.a and xB_val <= 0.5 * self.a:
                alpha = (self.a + 2*xA_val*sign)/(xB_val - xA_val)*sign
                r1 = np.array([(xB_val - xA_val)/3-0.5*self.a*sign, 0, -self.c/6])
                r2 = np.array([(xB_val-0.5*self.a*sign)/2, 0, 0])
                cg = 1 / (1+alpha) * r1 + alpha/(1+alpha) * r2
                # print(cg)
                cg[0] = cg[0] + self.origin[0]
                cg[1] = cg[1] + self.origin[1]
                cg[2] = cg[2] + self.origin[2]
                return cg
        
        # case 3
        zA = sp.symbols('zA')
        eqns = [zA- zB- self.a*np.tan(self.theta)]
        eqns.append(self.mass/self.M-0.5-(zA+zB)/(2*self.c))
        S = sp.solve(eqns, [zA, zB], dict=True, real=True)
        if len(S) == 1:
            # print('case 3')
            zA_val = S[0][zA]
            zB_val = S[0][zB]
            if zA_val >= -0.5 * self.c and zA_val <= 0.5 * self.c and zB_val >= -0.5 * self.c and zB_val <= 0.5 * self.c:
                alpha = (zA_val + 0.5 * self.c) / (zB_val + 0.5 * self.c)
                r1 = np.array([-0.5*self.a, 0, zA_val+zB_val-0.5*self.c])/3
                r2 = np.array([0.5*self.a, 0, zB_val - self.c])/3
                cg = alpha / (1+alpha) * r1 + 1/(1+alpha) * r2
                # print(cg)
                cg[0] = cg[0] + self.origin[0]
                cg[1] = cg[1] + self.origin[1]
                cg[2] = cg[2] + self.origin[2]
                return cg

        # case 4
        eqns = [np.tan(self.theta)*(xA+0.5*self.a*sign) - 0.5*self.c - zB]
        eqns = [self.mass/self.M - 0.5 * (0.5*self.a + xA*sign) * (zB + 0.5*self.c)/(self.a * self.c)]
        S = sp.solve(eqns, [xA, zB], dict=True, real=True)
        if len(S) == 1:
            # print('case 4')
            xA_val = S[0][xA]
            zB_val = S[0][zB]
            if xA_val >= -0.5 * self.a and xA_val <= 0.5 * self.a and zB_val >= -0.5 * self.c and zB_val <= 0.5 * self.c:
                cg = np.array([xA_val - self.a * sign, 0, zB_val- self.c])/3
                # print(cg)
                cg[0] = cg[0] + self.origin[0]
                cg[1] = cg[1] + self.origin[1]
                cg[2] = cg[2] + self.origin[2]
                return cg

        raise ValueError("No solution found when compute CG!")



if __name__ == '__main__':
    file_name = './data/oil_tank_params.json'
    fuel_feed_system = FuelFeedSystem(file_name=file_name)

    print(fuel_feed_system.topology_matrix)
    print(fuel_feed_system.tanks_directly_to_engine)
    print(fuel_feed_system.center_of_gravity())
    print(fuel_feed_system.output_fuel_mass_upper_bound)

    fuel_feed = np.array([1, 1, 1, 1, 1, 1])
    fuel_feed_system.step(fuel_feed)