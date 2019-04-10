'''
Example two qubit gate (sqwt swap based cphase)
'''

from qdot_wrappers import double_dot_hamiltonian
import numpy as np
import matplotlib.pyplot as plt

solver = double_dot_hamiltonian(2e9,2e9,1e12,1e12,0e1)


solver.awg_pulse_tc(4.3e9,1e-9,5e-9)
solver.awg_pulse_B_field(-10e6, 10e6, 10e-9, 35e-9)
solver.awg_pulse_tc(4.3e9,40e-9,44e-9)

dm = np.zeros([6,6], dtype=np.complex)
dm[1,1] = 1

solver.calc_time_evolution(dm, 0, 50e-9,25000)

print(solver.get_unitary())

solver.plot_pop()
solver.plot_expect()
solver.plot_ST_bloch_sphere()
plt.show()