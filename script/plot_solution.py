import numpy as np
import matplotlib.pyplot as plt
from coordinate_transformation import xyz2blh

def plot_trajectory(solution):
    """
    Plot the trajectory based on the provided solution data.
    :param solution: 2D array where each row contains [time, x, y, z] coordinates.
    """
    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    pos = np.zeros((nsol, 3))
    pos0 = np.zeros((nsol, 3))
    n = 0

    for i in range(nsol):
        # Extract positions (assuming the data format is [time, x, y, z])
        pos[n, :] = solution[i, 2:5]
        n += 1

    # Transform from XYZ to BLH (assuming xyz2blh is defined)
    _, Cne = xyz2blh(solution[0, 2:5])  # Transformation matrix at the first position

    j = 0
    for i in range(n):
        if np.dot(pos[i, 0:3], pos[i, 0:3]) <= 0:
            continue
        pos0[j, :] = np.dot(Cne, pos[i, 0:3])
        j += 1

    if j == 0:
        raise ValueError('Solution is empty!!!')

    if j < nsol:
        pos0 = pos0[:j, :]

    # Plot
    H = plt.gcf().get_size_inches() * plt.gcf().dpi  # Get screen size
    w, h = 600, 450
    x = (H[0] - w) / 2
    y = (H[1] - h) / 2

    # plt.plot(pos0[:, 0] - pos0[0, 0], pos0[:, 1] - pos0[0, 1], ':', linewidth=0.1, color=[0.5, 0.5, 0.5])
    plt.plot(pos0[:, 0] - pos0[0, 0], pos0[:, 1] - pos0[0, 1], '.b', markersize=3)

    plt.grid(True, linestyle='--', linewidth=1.0, color='k', alpha=0.3)
    plt.xlabel('E [m]', fontsize=12, family='Times New Roman')
    plt.ylabel('N [m]', fontsize=12, family='Times New Roman')
    plt.title('Trajectory', fontsize=12, family='Times New Roman')
    plt.axis('equal')
    plt.gca().ticklabel_format(style='plain', axis='x')  # Disable scientific notation for x-axis

    # 调用plt.show()时会阻塞程序执行，直到用户关闭当前图形窗口后才继续执行后续代码
    # plt.show()
    # 返回图形对象
    return plt.gcf()

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

    # 返回图形对象
    return fig

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
    fig, axes = plt.subplots(3, 1)
    Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

    # Plot East (E) velocity change
    axes[0].plot(time, vel0[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[0].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[0].set_ylabel('E [m/s]', fontsize=12, family='Times New Roman')
    axes[0].set_title('Velocity', fontsize=12, family='Times New Roman')

    # Plot North (N) velocity change
    axes[1].plot(time, vel0[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[1].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[1].set_ylabel('N [m/s]', fontsize=12, family='Times New Roman')

    # Plot Up (U) velocity change
    axes[2].plot(time, vel0[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[2].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[2].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    axes[2].set_ylabel('U [m/s]', fontsize=12, family='Times New Roman')

    # Remove scientific notation for axis labels
    for ax in axes:
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')

    # Adjust spacing
    plt.tight_layout()

    # 返回图形对象
    return fig

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
    fig, axes = plt.subplots(3, 1)
    Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

    # Plot East (E) velocity change
    axes[0].plot(time, att[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[0].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[0].set_ylabel('Pitch [deg]', fontsize=12, family='Times New Roman')
    axes[0].set_title('Attitude', fontsize=12, family='Times New Roman')

    # Plot North (N) velocity change
    axes[1].plot(time, att[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[1].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[1].set_ylabel('Roll [deg]', fontsize=12, family='Times New Roman')

    # Plot Up (U) velocity change
    axes[2].plot(time, att[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    axes[2].grid(True, linestyle='--', color='k', alpha=0.3)
    axes[2].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    axes[2].set_ylabel('Yaw [deg]', fontsize=12, family='Times New Roman')

    # Remove scientific notation for axis labels
    for ax in axes:
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')

    # Adjust spacing
    plt.tight_layout()

    # 返回图形对象
    return fig

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
