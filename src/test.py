import os
import sys
 
sys.path = sys.path[3:] + sys.path[:3]

import numpy as np
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool as ThreadPool
import time
import matlab.engine

# time_0 = time.time()
# cg_of_tanks_1 = []
# for i in range(1, 10):
#     time.sleep(1)
#     # cg_of_tanks_1.append(list(np.array([i,i+1,i+2])))
# time_1 = time.time()
# # print(np.array(cg_of_tanks_1))

# print ( "CPU的核数为: {}". format (cpu_count()))


# def process(i):
#     # cg_of_tanks_2.append(list(np.array([i,i+1,i+2])))
#     time.sleep(1)

# time_2 = time.time()
# cg_of_tanks_2 = []
# tanks_id = range(1, 10)
# pool = ThreadPool()
# pool.map(process, tanks_id)
# pool.close()
# pool.join()
# time_3 = time.time()

# # print(np.array(cg_of_tanks_2))

# print('time:', time_1 - time_0)
# print('parallel time:', time_3 - time_2)


# matlab call
# function c = CoM(mass, theta, shapeSize, oil_density, origin)

eng = matlab.engine.start_matlab()
d = eng.CoM(1000, 0, [1, 2, 3], 850, [0, 0, 0])
print(d)
