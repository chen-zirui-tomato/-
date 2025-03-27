import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 读取数据 
Z = np.loadtxt('../opt/data.txt')
ny, nx = Z.shape
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
X, Y = np.meshgrid(x, y)

# 将 Z=0 的点标记为 NaN（跳过不绘制）
Z_masked = np.where(Z == 0, np.nan, Z)  # 核心修改

# 绘制三维图
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 创建颜色映射（基于有效值范围）
norm = plt.Normalize(vmin=np.nanmin(Z_masked), vmax=np.nanmax(Z_masked))  # 忽略 NaN 计算范围
colors = plt.cm.rainbow(norm(Z_masked))  # 自动处理 NaN

# 绘制表面（会自动跳过 NaN 点）
ax.plot_surface(
    X, Y, Z_masked,
    rstride=1,
    cstride=1,
    facecolors=colors,
    shade=False  # 关闭光照避免 NaN 干扰
)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.savefig('../opt/plot.png')
plt.show(block=False)
input("Press Enter to continue...")
plt.close()