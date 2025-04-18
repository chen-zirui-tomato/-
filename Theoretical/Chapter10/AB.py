import numpy as np
import matplotlib.pyplot as plt

orders = [1, 2, 3, 4, 5]
theta = np.linspace(0, 2 * np.pi, 1000)
z = np.exp(1j * theta)

plt.figure(figsize=(10, 6))
plt.title("Adams-Bashforth 绝对稳定域 (p=1 到 5)")
plt.xlabel("Re(κ)")
plt.ylabel("Im(κ)")
plt.grid(True)

for p in orders:
    if p == 1:
        rho_coeff = [1, -1]        
        sigma_coeff = [0, 1]       
    elif p == 2:
        rho_coeff = [1, -1, 0]     
        sigma_coeff = [0, 3/2, -1/2]
    elif p == 3:
        rho_coeff = [1, -1, 0, 0]
        sigma_coeff = [0, 23/12, -16/12, 5/12]
    elif p == 4:
        rho_coeff = [1, -1, 0, 0, 0]
        sigma_coeff = [0, 55/24, -59/24, 37/24, -9/24]
    elif p == 5:
        rho_coeff = [1, -1, 0, 0, 0, 0]
        sigma_coeff = [0, 1901/720, -2774/720, 2616/720, -1274/720, 251/720]

    rho = np.polyval(rho_coeff, z)
    sigma = np.polyval(sigma_coeff, z)
    gamma = rho / sigma

    plt.plot(np.real(gamma), np.imag(gamma), label=f"p={p}")

plt.legend()
plt.show()