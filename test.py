import numpy as np
import time
import integral_lattice as il

strat_time = time.time()

t = 0
A = il.intergral_lattice(np.array([[2,t],[t,4]]))
#A = il.intergral_lattice(np.array([[3,1,1],[1,7,t],[1,t,13]]))
#A = il.intergral_lattice(np.array([[3,1,0],[1,3,t],[0,t,4]]))
#A = il.intergral_lattice(np.array([[3,1,0,0],[1,7,0,0],[0,0,2,0],[0,0,0,2]]))
#A = il.intergral_lattice(np.array([[2,0,0,0,0],[0,2,0,0,0],[0,0,2,0,0],[0,0,0,2,0],[0,0,0,0,1]]))


print(A.intersection_form)
print(A.disc)
A.dual_group_generator_frac()



duration = round(time.time() - strat_time,2)
print("Time cost:{}".format(duration))