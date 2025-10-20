import numpy as np
import matplotlib.pyplot as plt
from xyz2blh import xyz2blh

def plot_velocity(solution):
    """
    Plot the velocity changes (E, N, U) over GPS time.
    :param solution: A numpy array or list where each row contains [time, x, y, z].
    """
    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    vel = np.zeros((nsol, 3))
    vel0 = np.zeros((nsol, 3))
    time = np.zeros(nsol)
    n = 0

    for i in range(nsol):
        time[n] = solution[i, 1] #GPS time
        vel[n, :] = solution[i, 6:9]  # Assuming solution[i] has [time, x, y, z]
        n += 1

    # Convert from ECEF to ENU frame (Cne is the transformation matrix)
    _, Cne = xyz2blh(solution[0, 2:5])  # Assuming xyz2blh is defined elsewhere

    j = 0
    for i in range(n):
        if np.dot(vel[i, 0:3], vel[i, 0:3]) <= 0:
            continue
        vel0[j, :] = np.dot(Cne, vel[i, 0:3])
        j += 1

    if j == 0:
        raise ValueError('Solution is empty!!!')

    if j < nsol:
        time = time[:j]
        vel0 = vel0[:j, :]

    # Set up the figure window
    fig, axes = plt.subplots(3, 1, figsize=(10, 6))
    scrsz = plt.gcf().get_size_inches() * plt.gcf().dpi  # Screen size (Width x Height)
    Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

    # Plot East (E) velocity change
    axes[0].plot(time, vel0[:, 0], color=Fcolor[0], linestyle='-', marker='.')
    axes[0].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[0].set_ylabel('E [m/s]', fontsize=12, family='Times New Roman')
    axes[0].tick_params(axis='y', labelsize=12)
    axes[0].tick_params(axis='x', labelsize=12)
    axes[0].set_title('Velocity', fontsize=12, family='Times New Roman')

    # Plot North (N) velocity change
    axes[1].plot(time, vel0[:, 1], color=Fcolor[1], linestyle='-', marker='.')
    axes[1].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[1].set_ylabel('N [m/s]', fontsize=12, family='Times New Roman')
    axes[1].tick_params(axis='y', labelsize=12)
    axes[1].tick_params(axis='x', labelsize=12)

    # Plot Up (U) velocity change
    axes[2].plot(time, vel0[:, 2], color=Fcolor[2], linestyle='-', marker='.')
    axes[2].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[2].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    axes[2].set_ylabel('U [m/s]', fontsize=12, family='Times New Roman')
    axes[2].tick_params(axis='y', labelsize=12)
    axes[2].tick_params(axis='x', labelsize=12)

    # Remove scientific notation for axis labels
    for ax in axes:
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')

    # Adjust spacing
    plt.tight_layout()

    # Show the plot
    plt.show()
