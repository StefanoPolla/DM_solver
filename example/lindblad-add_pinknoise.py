import c_solver.ME_solver as me
import numpy as np
import matplotlib.pyplot as plt
import time 

samples = int(1e2)
evo_steps = int(1e3)
white_ampl = 0.5
pink_ampl = 0.5
total_t = 10

Z = np.array([[1,0],[0,-1]], dtype=complex)
X = np.array([[0,1],[1, 0]], dtype=complex)
Z0 = np.array([[0,0],[0,1]], dtype=complex)
Z1 = Z+Z0
pop = (np.array([Z0, Z1]), ['0','1'])

# Lindblad calculation and plotting
lindstart = time.time()

s = me.Lindblad(2) # this is a copy of VonNeumann with added support for lindblad operators
s.add_H0(Z) # works exactly like for VonNeumann
s.add_lindbladian(X, white_ampl**2) # produces the same results of a white noise, but takes as second input the noise variance instead of the amplitude (var = ampl^2)
pink = me.noise_py()
pink.init_pink(X, pink_ampl, 1)
s.add_noise_obj(pink)
s.set_number_of_evalutions(samples)

lindmid = time.time()

s.calculate_evolution(Z0, 0, total_t, evo_steps)

lindend = time.time()

s.plot_expectation(*pop, 1)

# VonNeumann calculation and plotting

vnstart = time.time()

v = me.VonNeumann(2)
v.add_H0(Z)
white = me.noise_py()
white.init_white(X, white_ampl)
v.add_noise_obj(white)
pink = me.noise_py()
pink.init_pink(X, pink_ampl, 1)
v.add_noise_obj(pink)
v.set_number_of_evalutions(samples)

vnmid = time.time() 

v.calculate_evolution(Z0, 0, total_t, evo_steps)

vnend = time.time()

v.plot_expectation(*pop, 2)

print('\n Measuring the time of the script:')
print(f'Lindblad approach   -- total time = {(lindend-lindstart):.5}s of which {(lindend-lindmid):.5}s of core operations.')
print(f'VonNeumann approach -- total time = {(vnend-vnstart):.5}s of which {(vnend-vnmid):.5}s of core operations.')
print(f'\n Both approaches were run on {samples} Monte Carlo samples,'+
	  f'with respective average evolution run time per sample {(lindend-lindmid)/samples:.5}ms and {(vnend-vnmid)/samples:.5}ms.')
print(f'In both cases, the number of evolution steps are {evo_steps}.')

plt.show()
