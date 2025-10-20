import numpy as np
import matplotlib.pyplot as plt

def plot_attitude(solution):
    """
    Plot the attitude over GPS time.
    :param solution: A numpy array or list where each row contains [time, x, y, z].
    """
    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    att = np.zeros((nsol, 3))
    time = np.zeros(nsol)
    n = 0

    for i in range(nsol):
        time[n] = solution[i, 1] #GPS time
        att[n, :] = solution[i, 9:12]  # Assuming solution[i] has [time, x, y, z]
        n += 1

    # Set up the figure window
    fig, axes = plt.subplots(3, 1, figsize=(10, 6))
    scrsz = plt.gcf().get_size_inches() * plt.gcf().dpi  # Screen size (Width x Height)
    Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

    # Plot East (E) velocity change
    axes[0].plot(time, att[:, 0], color=Fcolor[0], linestyle='-', marker='.')
    axes[0].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[0].set_ylabel('Pitch [deg]', fontsize=12, family='Times New Roman')
    axes[0].tick_params(axis='y', labelsize=12)
    axes[0].tick_params(axis='x', labelsize=12)
    axes[0].set_title('Attitude', fontsize=12, family='Times New Roman')

    # Plot North (N) velocity change
    axes[1].plot(time, att[:, 1], color=Fcolor[1], linestyle='-', marker='.')
    axes[1].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[1].set_ylabel('Roll [deg]', fontsize=12, family='Times New Roman')
    axes[1].tick_params(axis='y', labelsize=12)
    axes[1].tick_params(axis='x', labelsize=12)

    # Plot Up (U) velocity change
    axes[2].plot(time, att[:, 2], color=Fcolor[2], linestyle='-', marker='.')
    axes[2].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[2].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    axes[2].set_ylabel('Yaw [deg]', fontsize=12, family='Times New Roman')
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
