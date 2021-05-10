import matplotlib.pyplot as plt
import numpy as np

# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D

plt.ion()

box_size = float(input())

fig = plt.figure()
ax = plt.subplot(121, projection='3d', aspect='equal')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim3d(0, box_size)
ax.set_ylim3d(0, box_size)
ax.set_zlim3d(0, box_size)
# The values of s seems good enough for our case
# (found by hit and trial)
if box_size == 10.0:
    scat = ax.scatter([], [], [], s=(500))
elif box_size == 20.0:
    scat = ax.scatter([], [], [], s=(100))
elif box_size == 50.0:
    scat = ax.scatter([], [], [], s=(10))
elif box_size == 100.0:
    scat = ax.scatter([], [], [], s=(1))
plt.grid()

ax2 = plt.subplot(122)

# Maximize window, works only on Windows OS
plt.get_current_fig_manager().window.state('zoomed')

prev_nPoints = -1

while True:
    xs = []
    ys = []
    zs = []
    colors = []

    try:
        nPoints = int(input())
        # This script blocks output of the C++ program
        # import time
        # print('zzz')
        # time.sleep(5000)
    except EOFError:
        plt.show(block=True)
        fig.savefig("lastrun.png")
        exit()

    if nPoints != prev_nPoints:
        plt.cla()
        # ax2.set_xlabel('Monte Carlo Steps')
        ax2.set_xlabel('Number of trials')
        ax2.set_ylabel('Energy (in units of epsilon)')
        plt.grid()

    prev_nPoints = nPoints

    # mcs, energy = list(map(float, input().split()))
    trial_num, energy = list(map(float, input().split()))

    for _ in range(nPoints):
        inp = list(map(float, input().split()))
        xs.append(inp[0])
        ys.append(inp[1])
        zs.append(inp[2])

    colors = list(
        map(lambda x: int(x) / (nPoints-(nPoints != 1)), input().split()))

    if plt.get_fignums():
        # expects marker size in points**2.
        # formula: (2*rpix/fig.dpi*72)**2
        # Dpix = (ax.transData.transform([1.0,1.0]) - ax.transData.transform((0.0,0.0)))
        # print(Dpix)
        # Dpix = Dpix[0]
        # print(fig.dpi)
        # print(((Dpix/fig.dpi)*72))
        # ax.scatter(xs, ys, zs, c='r', s=(Dpix/fig.dpi)**2)

        # Temporary fix
        #ax.scatter(xs, ys, zs, c='r', s=(300))
        scat._offsets3d = (xs, ys, zs)
        # Think about Python linter issue on this line
        scat._facecolor3d = (list(map(plt.cm.gist_ncar, colors)))
        # plt.plot(mcs, energy, 'bo')
        plt.plot(trial_num, energy, 'bo')
        plt.pause(0.00001)
    else:
        fig.savefig("lastrun.png")
        exit()
