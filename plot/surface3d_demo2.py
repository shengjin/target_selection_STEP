from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = 20 * np.outer(np.cos(u), np.sin(v))
y = 20 * np.outer(np.sin(u), np.sin(v))
z = 20 * np.outer(np.ones(np.size(u)), np.cos(v))

ax.plot_surface(x, y, z, alpha=0.1, rstride=9, cstride=9, color='w')


target_array = np.genfromtxt('tgplot.csv', delimiter=",")
ax.scatter(target_array[:,18],target_array[:,19],target_array[:,20],",")

plt.show()

