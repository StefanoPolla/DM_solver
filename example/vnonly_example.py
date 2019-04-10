import c_solver.ME_solver as me
import numpy as np
import matplotlib.pyplot as plt
import time 

vn_samples = int(1e1)
evo_steps = int(2e4)
noise_ampl = 0.5
total_t = 10

Z = np.array([[1,0],[0,-1]], dtype=complex)
X = np.array([[0,1],[1, 0]], dtype=complex)
Z0 = np.array([[0,0],[0,1]], dtype=complex)
Z1 = Z+Z0
pop = (np.array([Z0, Z1]), ['0','1'])

# VonNeumann calculation and plotting

vnstart = time.time()

v = me.VonNeumann(2)
v.add_H0(Z)
white = me.noise_py()
white.init_white(X, noise_ampl)
v.add_noise_obj(white)
v.set_number_of_evalutions(vn_samples)

vnmid = time.time() 

v.calculate_evolution(Z0, 0, total_t, evo_steps)

vnend = time.time()

v.plot_expectation(*pop, 2)


print(f'VonNeumann approach -- total time = {(vnend-vnstart):.5}s of which {(vnend-vnmid)*1e3:.2}ms for evolution.')
print(f'\n The VonNeumann was run on {vn_samples} Monte Carlo samples, with average evolution run time per sample {(vnend-vnmid)/vn_samples*1e3:.2}ms.')
print(f'the number of evolution steps are {evo_steps}.')

plt.show()
