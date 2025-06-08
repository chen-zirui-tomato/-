import numpy as np
import matplotlib.pyplot as plt

# Parameters
h = 0.05
x = np.arange(-2, 25 + h, h)
T = 17
a = 1.0
mu_values = [0.8, 1.0]
k_values = [mu * h for mu in mu_values]
initial = np.exp(-20 * (x - 2)**2) + np.exp(-(x - 5)**2)

# Define numerical methods
def lax_friedrichs(u0, a, h, k, T):
    u = u0.copy()
    steps = int(T / k)
    mu = a * k / h
    for _ in range(steps):
        u[1:-1] = 0.5 * (u[2:] + u[:-2]) - mu/2 * (u[2:] - u[:-2])
    return u

def upwind(u0, a, h, k, T):
    u = u0.copy()
    steps = int(T / k)
    mu = a * k / h
    for _ in range(steps):
        u[1:] = u[1:] - mu * (u[1:] - u[:-1])
    return u

def lax_wendroff(u0, a, h, k, T):
    u = u0.copy()
    steps = int(T / k)
    mu = a * k / h
    for _ in range(steps):
        u[1:-1] = u[1:-1] - 0.5 * mu * (u[2:] - u[:-2]) + 0.5 * mu**2 * (u[2:] - 2*u[1:-1] + u[:-2])
    return u

def leapfrog(u0, a, h, k, T):
    u_prev = u0.copy()
    steps = int(T / k)
    mu = a * k / h
    u = lax_friedrichs(u0.copy(), a, h, k, k)  # 用 Lax-Friedrichs 初始化
    for _ in range(1, steps):
        u_new = u_prev.copy()
        u_new[1:-1] = u_prev[1:-1] - mu * (u[2:] - u[:-2])
        u_prev, u = u, u_new
    return u

def beam_warming(u0, a, h, k, T):
    u = u0.copy()
    steps = int(T / k)
    mu = a * k / h
    for _ in range(steps):
        u[2:] = u[2:] - 0.5 * mu * (3*u[2:] - 4*u[1:-1] + u[:-2]) + 0.5 * mu**2 * (u[2:] - 2*u[1:-1] + u[:-2])
    return u

# map
results = {}
schemes = {
    'leapfrog': leapfrog,
    'upwind': upwind,
    'beam_warming': beam_warming,
    'lax_friedrichs': lax_friedrichs,
    'lax_wendroff': lax_wendroff,
}

for mu in mu_values:
    k = mu * h
    for name, method in schemes.items():
        results[(name, mu)] = method(initial.copy(), a, h, k, T)

# exact solution    
def exact_solution(x, T):
    return np.exp(-20 * ((x - T) - 2)**2) + np.exp(-((x - T) - 5)**2)

u_exact = exact_solution(x, T)

# Plot
fig, axs = plt.subplots(4, 2, figsize=(10, 12))
axs = axs.flatten()

axs[0].plot(x, initial, 'k')
axs[0].set_title("initial condition")
axs[0].grid(True)

plot_order = [
    ('leapfrog', 0.8), ('upwind', 0.8), ('beam_warming', 0.8),
    ('lax_friedrichs', 0.8), ('lax_wendroff', 0.8),
    ('lax_wendroff', 1.0), ('leapfrog', 1.0)
]

for i, (name, mu) in enumerate(plot_order, start=1):
    axs[i].plot(x, u_exact, 'k-', linewidth=1)
    axs[i].plot(x, results[(name, mu)], 'b:', marker='o', markersize=0.5)
    axs[i].set_title(f"{name.replace('_', '-')} (k = {mu}h)")
    axs[i].grid(True)

plt.tight_layout()
plt.show()
