import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("../test/output.txt")
x = data[:, 0]
y = data[:, 1]

plt.figure(figsize=(6, 6))
plt.plot(x, y, label="approx")

# # 取每隔100个点一个作为控制点
# idx = np.arange(0, len(x), 100)
# plt.scatter(x[idx], y[idx], s=10, c='red', label="numerical answer")

# # 标注控制点编号
# for i in range(len(idx)):
#     plt.text(x[idx[i]], y[idx[i]], str(idx[i]), fontsize=8, color='blue')

plt.axis('equal')
# plt.gca().invert_yaxis()  # Y轴从上往下递增
plt.legend()
plt.title('IVPs')
plt.grid(True)
plt.show()