
import numpy as np
import matplotlib.pyplot as plt
import ttisect as T


# define triangles
# case -1
# t1 = np.asarray([[0,0,0], [0,1,0], [1,0,0]], dtype=np.float64)
# t2 = np.asarray([[1,0,0], [0,1,0], [1,1,0]], dtype=np.float64)

# case - 2
# t1 = np.asarray([[0,0,0], [0,1,0], [1,0,0]], dtype=np.float64)
# t2 = np.asarray([[0,0,0], [0,1,0], [-1,0,0]], dtype=np.float64)

t1 = np.asarray([[0,0,0], [0,1,0], [1,0,0]], dtype=np.float64)
t2 = np.asarray([[2,0,0], [2,1,0], [1,0,0]], dtype=np.float64)

# t2 = t2 + 0.1
# plot the triangles
def plot_triangles(fig, ax, t1, clr):
    ax.plot([t1[0][0], t1[1][0]],[t1[0][1], t1[1][1]], 
            '-', color=clr)
    ax.plot([t1[0][0], t1[2][0]],[t1[0][1], t1[2][1]], 
            '-', color=clr)
    ax.plot([t1[1][0], t1[2][0]],[t1[1][1], t1[2][1]], 
            '-', color=clr)
    # ax.plot(t1[2][:-1], t1[0][:-1], '-', color=clr)

fig, ax = plt.subplots()
plot_triangles(fig, ax, t1, 'k')
plot_triangles(fig, ax, t2, 'b')
plt.show()

print(T.tri_tri_isect(t1[0],t1[2],t1[1],
    t2[1],t2[0],t2[2]))
