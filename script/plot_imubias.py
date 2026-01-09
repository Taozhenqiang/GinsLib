import numpy as np
import matplotlib.pyplot as plt

def plot_imubias(solution):
    """
    Plot gyroscope and accelerometer biases over GPS time.
    :param solution: A numpy array where each row contains [time, ..., bgx, bgy, bgz, bax, bay, baz].
    """
    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    bg = np.zeros((nsol, 3))  # Gyroscope biases (bgx, bgy, bgz)
    ba = np.zeros((nsol, 3))  # Accelerometer biases (bax, bay, baz)
    time = np.zeros(nsol)

    for n in range(nsol):
        time[n] = solution[n, 1]  # GPS time
        bg[n, :] = solution[n, 12:15]  # Gyroscope biases (deg/h)
        ba[n, :] = solution[n, 15:18]  # Accelerometer biases (ug) 

    # Plot gyroscope biases
    Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

    # Gyroscope bias plot
    fig_bg, ax_bg = plt.subplots(3, 1)
    
    ax_bg[0].plot(time, bg[:, 0], color=Fcolor[0],  linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax_bg[0].grid(True, linestyle='--', color='k', alpha=0.3)
    ax_bg[0].set_ylabel('bgx [deg/h]', fontsize=12, family='Times New Roman')
    ax_bg[0].set_title('Gyroscope bias', fontsize=12, family='Times New Roman')

    ax_bg[1].plot(time, bg[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax_bg[1].grid(True, linestyle='--', color='k', alpha=0.3)
    ax_bg[1].set_ylabel('bgy [deg/h]', fontsize=12, family='Times New Roman')

    ax_bg[2].plot(time, bg[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax_bg[2].grid(True, linestyle='--', color='k', alpha=0.3)
    ax_bg[2].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax_bg[2].set_ylabel('bgz [deg/h]', fontsize=12, family='Times New Roman')

    # Remove scientific notation for axis labels
    for ax in ax_bg:
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')

    plt.tight_layout()


    # Plot accelerometer biases
    fig_ba, ax_ba = plt.subplots(3, 1)

    # Accelerometer bias plot
    ax_ba[0].plot(time, ba[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax_ba[0].grid(True, linestyle='--', color='k', alpha=0.3)
    ax_ba[0].set_ylabel('bax [ug]', fontsize=12, family='Times New Roman')
    ax_ba[0].set_title('Accelerometer bias', fontsize=12, family='Times New Roman')

    ax_ba[1].plot(time, ba[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax_ba[1].grid(True, linestyle='--', color='k', alpha=0.3)
    ax_ba[1].set_ylabel('bay [ug]', fontsize=12, family='Times New Roman')

    ax_ba[2].plot(time, ba[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax_ba[2].grid(True, linestyle='--', color='k', alpha=0.3)
    ax_ba[2].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax_ba[2].set_ylabel('baz [ug]', fontsize=12, family='Times New Roman')

    # Remove scientific notation for axis labels
    for ax in ax_ba:
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')

    plt.tight_layout()
    
    # 返回两个图形对象
    return fig_bg, fig_ba