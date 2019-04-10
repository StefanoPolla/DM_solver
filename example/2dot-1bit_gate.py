'''
Example single qubit gate
'''

from qdot_wrappers import double_dot_hamiltonian
import numpy as np
import matplotlib.pyplot as plt

test_single_qubit_gate = double_dot_hamiltonian(2e9,2.05e9,2e12,2.1e12,1e9)

# test_single_qubit_gate.mw_pulse(1.7e9,0,2e6,0e-9,500e-9, RWA= True)
# test_single_qubit_gate.mw_pulse(2.15e9,0,5e6,0e-9,500e-9, RWA= True)
test_single_qubit_gate.mw_pulse(2.05e9, 0, 2e6, 0e-9, 500e-9)


DM = np.zeros([6,6], dtype=np.complex)
DM[1,1] = 1

test_single_qubit_gate.calc_time_evolution(DM, 0, 600e-9, 10000)

test_single_qubit_gate.plot_pop()
plt.show()
