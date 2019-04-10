from qdot_wrappers import two_qubit_sim
import numpy as np
import matplotlib.pyplot as plt



''' Decreasing t_gate increases offresonant component (crosstalk) '''
t_gate = 10e-9 # good results with 100ns, at 10ns crosstalk is strong
bloch_resolution = 10000# main running time limitation

B1 = 2e9
B2 = 2.100e9 
mw_f = B1
rabi_f = 1 / (2 * t_gate) # amplitude tuned to half a rabi oscillation

pulse_start = t_gate * 0.1
tot_time = t_gate + 2*pulse_start

sim = two_qubit_sim(B1, B2)
sim.mw_pulse(mw_f, 0, rabi_f, pulse_start, pulse_start + t_gate)

DM = np.zeros([4,4], dtype=np.complex)
DM[0,0] = 1
sim.calc_time_evolution(DM, 0, tot_time, 40000)

sim.plot_pop()
# qutip.Bloch.show() blocks the program. choose which bloch sphere to see:
# first bit - the one that is supposed to flip
# sim.plot_qubit_bloch_sphere(bloch_resolution, 0, fig=plt.figure(2))
# second bit - the one shifted by crosstalk
sim.plot_qubit_bloch_sphere(bloch_resolution, 1)