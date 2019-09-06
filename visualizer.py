import matplotlib.pyplot as plt

# colormap = {
#     'empty': 'lightgrey',
#     'H'    : 'yellow',
#     'CO'   : 'brown',
#     'C'    : 'black',
#     'CH'   : 'lightgreen',
#     'CH2'  : 'green',
#     'CH3'  : 'darkgreen'
# }

colormap = {
    'empty': 'lightgrey',
    'H'    : 'yellow',
    'CO'   : 'brown',
    'C'    : 'black',
    'CH'   : 'cyan',
    'CH2'  : 'blue',
    'CH3'  : 'darkblue'
}

plt.ion()

lattice_size = int(input())

# Initialize the lattice
lattice = ['empty'] * (lattice_size * lattice_size)

fig = plt.figure()

# Add first plot (for lattice)
ax = plt.subplot(121, aspect='equal')

# Set axes' limits
ax.axis([-0.5, lattice_size, -0.5, (lattice_size-1)*(3**0.5/2) + 0.5])

# Invert the y-axis
plt.gca().invert_yaxis()

# Hide the axes
ax.axis('off')

# Initialize plot data
xs = [i % lattice_size + 0.5 * ((i//lattice_size) & 1)
      for i in range(lattice_size*lattice_size)]
ys = [(i//lattice_size) * (3**0.5/2) for i in range(lattice_size*lattice_size)]

# Display the legend
for _label, _color in colormap.items():
    ax.scatter([], [], color=_color, label=_label)
ax.legend()

# Plot the lattice on screen, with default color (to be changed later)
# The value of s seems good enough for our case
# (found by hit and trial)
# For 10x10, s=1400 and for 100x100, s=10 are found to work best
# scat = ax.scatter(xs, ys, s=1400) # For lattice size 10x10
scat = ax.scatter(xs, ys, s=10) # For lattice size 100x100

# Add second plot (for plotting dependant variables)
ax2 = plt.subplot(122)
ax2.grid()

# Maximize window, works only on Windows OS
plt.get_current_fig_manager().window.state('zoomed')

# Initialize the colors list
colors = [''] * (lattice_size * lattice_size)

prev_mole_fraction = -1

while True:
    try:
        mole_fraction = float(input())
        mcs = float(input())
        if mole_fraction != prev_mole_fraction:
            plt.cla()
            ax2.set_xlabel('Monte Carlo Steps')
            ax2.set_ylabel(f'Hydrocarbon production for mole fraction {mole_fraction}')
            ax2.grid()
        prev_mole_fraction = mole_fraction
        # This script blocks output of the C++ program
        # import time
        # print('zzz')
        # time.sleep(5000)
    except EOFError:
        plt.show(block=True)
        fig.savefig("lastrun.png")
        exit()

    lattice = list(input().split())
    hydrocarbon_production = int(input())

    if plt.get_fignums():
        colors = [colormap[key] for key in lattice]
        scat.set_facecolors(colors)
        plt.plot(mcs, hydrocarbon_production, 'bo')
        plt.pause(0.00001)
    else:
        fig.savefig("lastrun.png")
        exit()
