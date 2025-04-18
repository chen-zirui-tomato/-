import numpy as np
import matplotlib.pyplot as plt

orders = [3, 4, 5]
theta = np.linspace(0, 2 * np.pi, 1000)
z = np.exp(1j * theta)

plt.figure(figsize=(10, 6))
plt.title("Adams-Moulton 绝对稳定域 (p=3 到 5)")
plt.xlabel("Re(κ)")
plt.ylabel("Im(κ)")
plt.grid(True)

for p in orders:
    if p == 3:
        rho_coeff = [1, -1, 0, 0]
        sigma_coeff = [5/12, 2/3, -1/12, 0]
    elif p == 4:
        rho_coeff = [1, -1, 0, 0, 0]
        sigma_coeff = [3/8, 19/24, -5/24, 1/24, 0]
    elif p == 5:
        rho_coeff = [1, -1, 0, 0, 0, 0]
        sigma_coeff = [251/720, 646/720, -264/720, 106/720, -19/720, 0]
    
    rho = np.polyval(rho_coeff, z)
    sigma = np.polyval(sigma_coeff, z)
    gamma = rho / sigma

    plt.plot(np.real(gamma), np.imag(gamma), label="p={}".format(p))

plt.legend()
plt.show()

