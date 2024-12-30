Λ = 3300
alpha =2.1e-9
β1 = 4e-7
β2 = 0.005
sigma1 = 0.0005
sigma2 = 0.1
sigma3 = 0.005
r1 = 0.04
r2 = 0.002
d1 = 0.000015
d2 = 0.0032

# Initialize variables with initial conditions
S = 75000000
E = 225
Q = 800
I = 58
R = 30000000
t=0

delta_t = 0.1  # Time step (day)
T = 60    # Total duration of the simulation

def ds_dt(S, E, Q, I, R):
    return Λ - alpha * S * E - sigma1 * S - β1 * S - d1 * S 

def de_dt(S, E, Q, I, R):
    return alpha * S * E - r1 * E - β2 * E - d1 * E

def dq_dt(S, E, Q, I, R):
    return β1 * S + β2 * E - r2 * Q - sigma2 * Q - d1 * Q

def di_dt(S, E, Q, I, R):
    return r1 * E + r2 * Q - sigma3 * I - d1 * I - d2 * I

def dr_dt(S, E, Q, I, R):
    return sigma1 * S + sigma2 * Q + sigma3 * I - d1 * R


def runge_kutta_step(S, E, Q, I, R, h):
    k1s = ds_dt(S, E, Q, I, R) * h
    k1e = de_dt(S, E, Q, I, R) * h
    k1q = dq_dt(S, E, Q, I, R) * h
    k1i = di_dt(S, E, Q, I, R) * h
    k1r = dr_dt(S, E, Q, I, R) * h

    k2s = ds_dt(S + 0.5 * k1s, E + 0.5 * k1e, Q + 0.5 * k1q, I + 0.5 * k1i, R + 0.5 * k1r) * h
    k2e = de_dt(S + 0.5 * k1s, E + 0.5 * k1e, Q + 0.5 * k1q, I + 0.5 * k1i, R + 0.5 * k1r) * h
    k2q = dq_dt(S + 0.5 * k1s, E + 0.5 * k1e, Q + 0.5 * k1q, I + 0.5 * k1i, R + 0.5 * k1r) * h
    k2i = di_dt(S + 0.5 * k1s, E + 0.5 * k1e, Q + 0.5 * k1q, I + 0.5 * k1i, R + 0.5 * k1r) * h
    k2r = dr_dt(S + 0.5 * k1s, E + 0.5 * k1e, Q + 0.5 * k1q, I + 0.5 * k1i, R + 0.5 * k1r) * h

    k3s = ds_dt(S + 0.5 * k2s, E + 0.5 * k2e, Q + 0.5 * k2q, I + 0.5 * k2i, R + 0.5 * k2r) * h
    k3e = de_dt(S + 0.5 * k2s, E + 0.5 * k2e, Q + 0.5 * k2q, I + 0.5 * k2i, R + 0.5 * k2r) * h
    k3q = dq_dt(S + 0.5 * k2s, E + 0.5 * k2e, Q + 0.5 * k2q, I + 0.5 * k2i, R + 0.5 * k2r) * h
    k3i = di_dt(S + 0.5 * k2s, E + 0.5 * k2e, Q + 0.5 * k2q, I + 0.5 * k2i, R + 0.5 * k2r) * h
    k3r = dr_dt(S + 0.5 * k2s, E + 0.5 * k2e, Q + 0.5 * k2q, I + 0.5 * k2i, R + 0.5 * k2r) * h

    k4s = ds_dt(S + k3s, E + k3e, Q + k3q, I + k3i, R + k3r) * h
    k4e = de_dt(S + k3s, E + k3e, Q + k3q, I + k3i, R + k3r) * h
    k4q = dq_dt(S + k3s, E + k3e, Q + k3q, I + k3i, R + k3r) * h
    k4i = di_dt(S + k3s, E + k3e, Q + k3q, I + k3i, R + k3r) * h
    k4r = dr_dt(S + k3s, E + k3e, Q + k3q, I + k3i, R + k3r) * h

    S += (k1s + 2 * k2s + 2 * k3s + k4s) / 6
    E += (k1e + 2 * k2e + 2 * k3e + k4e) / 6
    Q += (k1q + 2 * k2q + 2 * k3q + k4q) / 6
    I += (k1i + 2 * k2i + 2 * k3i + k4i) / 6
    R += (k1r + 2 * k2r + 2 * k3r + k4r) / 6

    return S, E, Q, I, R

import matplotlib.pyplot as plt

# Lists to store the values 
time_steps = []
susceptible_values = []
exposed_values = []
quarantined_values = []
infected_values = []
recovered_values = []

# Main simulation loop 
while t < T:
    time_steps.append(t)
    susceptible_values.append(S)
    exposed_values.append(E)
    quarantined_values.append(Q)
    infected_values.append(I)
    recovered_values.append(R)
    
    S, E, Q, I, R = runge_kutta_step(S, E, Q, I, R, delta_t)
    t += delta_t

# Plotting the graphs
"""plt.figure(figsize=(10, 6))
plt.plot(time_steps, susceptible_values, label='Susceptible', color='blue')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Maharashtra')
plt.legend()
plt.show()
plt.plot(time_steps, exposed_values, label='Exposed', color='orange')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Maharashtra')
plt.legend()
plt.show()
plt.plot(time_steps, quarantined_values, label='Quarantined', color='green')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Maharashtra')
plt.legend()
plt.show()"""

plt.plot(time_steps, recovered_values, label='Recovered', color='purple')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Maharashtra')
plt.legend()
plt.show()
"""
plt.plot(time_steps, infected_values, label='Infected', color='red')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Maharashtra')
plt.legend()
plt.show()"""


