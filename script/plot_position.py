import numpy as np
import matplotlib.pyplot as plt
from xyz2blh import xyz2blh

def plot_position(solution):
    """
    Plot the position changes (E, N, U) over GPS time.
    :param solution: A numpy array or list where each row contains [time, x, y, z].
    """
    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    pos = np.zeros((nsol, 3))
    pos0 = np.zeros((nsol, 3))
    time = np.zeros(nsol)
    n = 0

    for i in range(nsol):
        time[n] = solution[i, 1] #GPS time
        pos[n, :] = solution[i, 2:5]  # Assuming solution[i] has [time, x, y, z]
        n += 1

    # Convert from ECEF to ENU frame (Cne is the transformation matrix)
    _, Cne = xyz2blh(solution[0, 2:5])  # Assuming xyz2blh is defined elsewhere

    j = 0
    for i in range(n):
        if np.dot(pos[i, 0:3], pos[i, 0:3]) <= 0:
            continue
        pos0[j, :] = np.dot(Cne, pos[i, 0:3])
        j += 1

    if j == 0:
        raise ValueError('Solution is empty!!!')

    if j < nsol:
        time = time[:j]
        pos0 = pos0[:j, :]

    # Calculate delta position
    delta = pos0[:, 0:3] - np.tile(pos0[0, 0:3], (j, 1))

    # Set up the figure window
    fig, axes = plt.subplots(3, 1)
    Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

    # Plot East (E) position change
    axes[0].plot(time, delta[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[0].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[0].set_ylabel('E [m]', fontsize=12, family='Times New Roman')
    axes[0].set_title('Position', fontsize=12, family='Times New Roman')

    # Plot North (N) position change
    axes[1].plot(time, delta[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[1].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[1].set_ylabel('N [m]', fontsize=12, family='Times New Roman')


    # Plot Up (U) position change
    axes[2].plot(time, delta[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[2].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[2].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    axes[2].set_ylabel('U [m]', fontsize=12, family='Times New Roman')


    # Remove scientific notation for axis labels
    for ax in axes:
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')

    # Adjust spacing
    plt.tight_layout()

    # Show the plot
    plt.show()
