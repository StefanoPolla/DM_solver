from qdot_wrappers import two_qubit_sim
import numpy as np
import matplotlib.pyplot as plt

B1 = 2e9
B2 = 2.100e9 

J = 10e6
start = 50e-9
t_gate = 1/J
end = 2*start + t_gate

sim = two_qubit_sim(B1, B2)

sim.RF_exchange(J, abs(B1-B2), start, start + t_gate, 0)

sim.pulsed_exchange(3*J/2, start, start + t_gate)

DM = np.zeros([4,4], dtype=np.complex)
DM[2,2] = 1

sim.calc_time_evolution(DM, 0, end, 40000)

U = sim.get_unitary()

U_wanted = np.array([
    [1,0,0,0],
    [0,0,1,0],
    [0,1,0,0],
    [0,0,0,1]
    ])

print(f'fidelity: {sim.get_fidelity(U_wanted):.2%} \n')

np.set_printoptions(precision=2)

print('target : \n', U_wanted, '\n')
print(f'U : \n', U, '\n' )
print(f'U * exp(3j/4*np.pi) : \n', 
      U*np.exp(3j/4*np.pi), '\n')

sim.plot_pop()
sim.plot_expect()
sim.plot_ST_bloch_sphere(500)