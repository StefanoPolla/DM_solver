import c_solver.ME_solver as me
import numpy as np
import matplotlib.pyplot as plt
import time 

evo_steps = int(2e4)
noise_ampl = 0.5 # TODO: did not yet check how to convert to T1
total_t = 10

Z = np.array([[1,0],[0,-1]], dtype=complex)
X = np.array([[0,1],[1, 0]], dtype=complex)
Z0 = np.array([[0,0],[0,1]], dtype=complex)
Z1 = Z+Z0
jumpdown = np.array([[0,0],[1, 0]], dtype=complex)
pop = (np.array([Z0, Z1]), ['0','1'])

# Lindblad calculation and plotting
lindstart = time.time()

s = me.Lindblad(2) # this is a copy of VonNeumann with added support for lindblad operators
s.add_H0(Z) # works exactly like for VonNeumann
s.add_lindbladian(jumpdown, noise_ampl**2) # with jumpdown, produces T1 decay. 

lindmid = time.time()

s.calculate_evolution(Z1, 0, total_t, evo_steps)

lindend = time.time()

s.plot_expectation(*pop, 1)

print('\n Measuring the time of the script:')
print(f'total time = {(lindend-lindstart):.5}s of which {(lindend-lindmid):.5}s of core operations.')
print(f'The number of evolution steps are {evo_steps}.')

plt.show()
