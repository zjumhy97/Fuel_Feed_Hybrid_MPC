import numpy as np
import sympy as sp

from parameters import configs

__all__ = ['FuelFeedSystem']


class FuelFeedSystem:
    def __init__(self, aircraft_mass:float, num_of_tanks:int, topology_matrix:np.ndarray) -> None:
        self.aircraft_mass = aircraft_mass
        self.num_of_tanks = num_of_tanks
        self.topology_matrix = topology_matrix
        self.mass_of_oil_in_tanks = 0

    def calculate_CG_of_aircraft(self, mass_of_oil_in_tanks:np.ndarray, ) -> np.ndarray:
        pass


class OilTank:
    def __init__(self, shape_size:np.ndarray, origin:np.ndarray) -> None:
        """
        Input:
        shape_size: three-demension vector, the size of the tank
        origin: the original point of the tank's coordinate frame w.r.t flight coordinate frame
        """
        self.shape_size = shape_size
        self.origin = origin
        self.a = shape_size[0]
        self.b = shape_size[1]
        self.c = shape_size[2]
        self.V = self.a * self.b * self.c
        self.M = configs.oil_density * self.V
    
    def _calculate_CG_of_single_tank(self, mass:np.ndarray, theta:float=0.0) -> np.ndarray:
        """
        Compute the center of gravity (CG) of a tank
        
        Input:
        mass: the mass of oil in the tank
        theta: the pitch angle of the flight

        Output:
        cg: the center of gravity (CG) of a tank
        """
        sign = 0 if theta < 0 else 1

        # case 1
        xA, zB = sp.symbols('xA zB')
        eqns = [zB - 0.5 * self.c + np.tan(theta) * (0.5 * self.a * sign - xA)]
        eqns.append((self.M - mass) / self.M - 0.5 * (0.5 * self.a - sign * xA) * (0.5 * self.c - zB) / (self.a * self.c))
        S = sp.solve(eqns, [xA, zB], dict=True, real=True)
        if len(S) == 1:
            # print('case 1')
            xA_val = S[0][xA]
            zB_val = S[0][zB]
            if xA_val >= -0.5 * self.a and xA_val <= 0.5 * self.a and zB_val >= -0.5 * self.c and zB_val <= 0.5 * self.c:
                cg = -(self.M - mass) / mass * np.array([(xA_val + self.a * sign) / 3, 0, (zB_val + self.c) / 3])
                cg[0] = cg[0] + self.origin[0]
                cg[1] = cg[1] + self.origin[1]
                cg[2] = cg[2] + self.origin[2]
                return cg

        # case 2
        xB = sp.symbols('xB')
        eqns = [np.tan(theta) * (xA - xB) + self.c]
        eqns.append(mass / self.M - 0.5 - 0.5/self.a*(xA+xB)*sign)
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
                print(cg)
                cg[0] = cg[0] + self.origin[0]
                cg[1] = cg[1] + self.origin[1]
                cg[2] = cg[2] + self.origin[2]
                return cg
        
        # case 3
        zA = sp.symbols('zA')
        eqns = [zA- zB- self.a*np.tan(theta)]
        eqns.append(mass/self.M-0.5-(zA+zB)/(2*self.c))
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
                print(cg)
                cg[0] = cg[0] + self.origin[0]
                cg[1] = cg[1] + self.origin[1]
                cg[2] = cg[2] + self.origin[2]
                return cg

        # case 4
        eqns = [np.tan(theta)*(xA+0.5*self.a*sign) - 0.5*self.c - zB]
        eqns = [mass/self.M - 0.5 * (0.5*self.a + xA*sign) * (zB + 0.5*self.c)/(self.a * self.c)]
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
    oil_tank_param_1 = {
        "shape_size":np.array([1.5, 0.9, 0.3]),
        "origin":np.array([8.91304348, 1.20652174, 0.61669004]),
    }
    mass_1 = 0.3 * configs.oil_density
    oil_tank_1 = OilTank(**oil_tank_param_1)
    res = oil_tank_1._calculate_CG_of_single_tank(mass=mass_1, theta=0.0)
    print('res=' , res)