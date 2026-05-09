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
        vel[n, :] = solution[i, 5:8]  # Assuming solution[i] has [time, x, y, z]
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
        att[n, :] = solution[i, 8:11]  # Assuming solution[i] has [time, x, y, z]
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
        if len(solution[n]) >= 18:
            bg[n, :] = solution[n, 15:18]  # Gyroscope biases (deg/h)
            ba[n, :] = solution[n, 18:21]  # Accelerometer biases (ug) 
        else:
            bg[n, :] = np.zeros(3)
            ba[n, :] = np.zeros(3)

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

def plot_nsatdop(solution):
    """
    Plot the number of satellites and pdop over GPS time.
    :param solution: A numpy array where each row contains [time, ..., nsv, dop].
    """

    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    nsat = np.zeros((nsol, 1))  # Number of satellites in view (NSV)
    pdop = np.zeros((nsol, 1))  # pdop (DOP)
    time = np.zeros(nsol)

    for n in range(nsol):
        time[n] = solution[n, 1]  # GPS time
        if len(solution[n]) >= 15:
            nsat[n, :] = solution[n, 12]  # nsat
            pdop[n, :] = solution[n, 14]  # pdop
        else:
            nsat[n, :] = np.zeros(1)
            pdop[n, :] = np.zeros(1)

    # Plot number of satellites and pdop
    Fcolor = ["#1a0dcf", "#14a959", "#ff6666"]

    # Number of satellites plot
    fig, ax = plt.subplots(2, 1)

    ax[0].plot(time, nsat, color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[0].grid(True, linestyle='--', color='k', alpha=0.3)
    ax[0].set_ylabel('NSAT', fontsize=12, family='Times New Roman')

    ax[1].plot(time, pdop, color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[1].grid(True, linestyle='--', color='k', alpha=0.3)
    ax[1].set_ylabel('PDOP', fontsize=12, family='Times New Roman')

    # Remove scientific notation for axis labels
    for ax in ax:
        ax.ticklabel_format(style='plain', axis='x')
        ax.ticklabel_format(style='plain', axis='y')

    plt.tight_layout()
    
    # 返回图形对象
    return fig

def plot_ratio(solution):
    """
    Plot the ratio over GPS time.
    :param solution: A numpy array where each row contains [time, ..., ratio].
    """

    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    ratio = np.zeros((nsol, 1))  # Ratio (ratio)
    time = np.zeros(nsol)  # GPS time  # ratio (ratio)

    for n in range(nsol):
        time[n] = solution[n, 1]  # GPS time
        if len(solution[n]) >= 13:
            ratio[n, :] = solution[n, 12]  # ratio
        else:
            ratio[n, :] = np.zeros(1)

    Fcolor = ["#37bdab"]

    # Plot ratio
    thres = 3.0  # Threshold for ratio
    ratio_idx = np.where(np.abs(ratio) > thres)[0]
    rate = len(ratio_idx)/nsol*100

    fig, ax = plt.subplots(1, 1)
    ax.plot(time, ratio, color=Fcolor[0], linestyle='none', linewidth=1.0, marker='.', markersize=2.5)
    ax.plot(time, thres*np.ones(nsol), color='r', linestyle='-', linewidth=1.0)
    ax.legend([f'Fix rate={rate:.2f}%(Thres={thres})'], fontsize=12)
    ax.grid(True, linestyle='--', color='k', alpha=0.3)
    ax.set_ylabel('Ratio', fontsize=12, family='Times New Roman')
    ax.set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax.ticklabel_format(style='plain', axis='x')
    ax.ticklabel_format(style='plain', axis='y')
    plt.tight_layout()  

    return fig

def plot_solflag(solution):
    """
    Plot the solution flag over GPS time.
    :param solution: A numpy array where each row contains [time, ..., solflag].
    """

    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')
    
    solflag = np.zeros((nsol, 1))  # Solution flag (solflag)
    time = np.zeros(nsol)  # GPS time  # solflag (solflag)
    
    for n in range(nsol):
        time[n] = solution[n, 1]  # GPS time
        if len(solution[n]) >= 12:
            solflag[n, :] = solution[n, 11]  # solflag
        else:
            solflag[n, :] = np.zeros(1)

    # extract solution flag    
    spp_idx = np.where(solflag == 5)[0]
    ppd_idx = np.where(solflag == 4)[0]
    float_idx = np.where(solflag == 2)[0]
    fix_idx = np.where(solflag == 1)[0]
    ins_idx = np.where(solflag == 7)[0]
    cons_idx = np.where(solflag == 8)[0]

    Fcolor = ["#14a959", "#ff6666", "#37bdb6", "#ff9900", "#0a4e21", "#1a0dcf",]

    # plot solution flag
    fig, ax = plt.subplots(1, 1)
    ax.plot(time[spp_idx], solflag[spp_idx], color=Fcolor[0], linestyle='none', marker='.', markersize=1.5)
    ax.plot(time[ppd_idx], solflag[ppd_idx], color=Fcolor[1], linestyle='none', marker='.', markersize=1.5)
    ax.plot(time[float_idx], solflag[float_idx], color=Fcolor[2], linestyle='none', marker='.', markersize=1.5)
    ax.plot(time[fix_idx], solflag[fix_idx], color=Fcolor[3], linestyle='none', marker='.', markersize=1.5)
    ax.plot(time[ins_idx], solflag[ins_idx], color=Fcolor[4], linestyle='none', marker='.', markersize=1.5)
    ax.plot(time[cons_idx], solflag[cons_idx], color=Fcolor[5], linestyle='none', marker='.', markersize=1.5)
    ax.legend([f'SPP={len(spp_idx)}({len(spp_idx)/nsol*100:.2f}%)', 
               f'PPD={len(ppd_idx)}({len(ppd_idx)/nsol*100:.2f}%)', 
               f'Float={len(float_idx)}({len(float_idx)/nsol*100:.2f}%)', 
               f'Fix={len(fix_idx)}({len(fix_idx)/nsol*100:.2f}%)',     
               f'INS={len(ins_idx)}({len(ins_idx)/nsol*100:.2f}%)', 
               f'Cons={len(cons_idx)}({len(cons_idx)/nsol*100:.2f}%)'], fontsize=12)
    ax.grid(True, linestyle='--', color='k', alpha=0.3)
    ax.set_ylabel('Solution Status', fontsize=12, family='Times New Roman')
    ax.set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax.ticklabel_format(style='plain', axis='x')
    ax.ticklabel_format(style='plain', axis='y')
    plt.tight_layout()  

    return fig

def plot_pvastd(solution):
    """
    Plot the position velocity standard deviation (PVASTD) over GPS time.
    :param solution: A numpy array where each row contains [time, ..., pvastd].
    """

    nsol = len(solution)
    if nsol == 0:
        raise ValueError('Solution is empty!!!')

    # Transform from XYZ to ENU frame
    _, Cne = xyz2blh(solution[0, 2:5])  # Transformation matrix at the first position

    posstd = np.zeros((nsol, 3))  # Position standard deviation (m) (posstd)
    velstd = np.zeros((nsol, 3))  # Velocity standard deviation (m/s) (velstd)
    attstd = np.zeros((nsol, 3))  # Attitude standard deviation (deg) (attstd)
    bgstd = np.zeros((nsol, 3))  # gyro bias standard deviation (deg/h) (bgstd)
    bastd = np.zeros((nsol, 3))  # accelerometer bias standard deviation (ug) (bastd)
    time = np.zeros(nsol)  # GPS time   

    for n in range(nsol):
        time[n] = solution[n, 1]  # GPS time
        if len(solution[n]) >= 15:
            pos_cov = Cne @ np.diagflat(solution[n, 21:24]**2) @ Cne.T # posstd
            posstd[n, :] = np.sqrt(np.diag(pos_cov))
            vel_cov = Cne @ np.diagflat(solution[n, 24:27]**2) @ Cne.T # velstd
            velstd[n, :] = np.sqrt(np.diag(vel_cov))
            attstd[n, :] = solution[n, 27:30]  # attstd
            bgstd[n, :] = solution[n, 30:33]  # bgstd
            bastd[n, :] = solution[n, 33:36]  # bastd
        else:
            posstd[n, :] = np.zeros(3)
            velstd[n, :] = np.zeros(3)
            attstd[n, :] = np.zeros(3)
            bgstd[n, :] = np.zeros(3)
            bastd[n, :] = np.zeros(3)

    # Plot pvt standard deviations
    Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

    # Gyroscope bias plot
    fig, ax = plt.subplots(3, 2)    

    ax[0, 0].plot(time, posstd[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[0, 0].plot(time, posstd[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[0, 0].plot(time, posstd[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[0, 0].grid(True, linestyle='--', color='k', alpha=0.3)
    ax[0, 0].legend(['E', 'N', 'U'], fontsize=12)
    ax[0, 0].set_ylabel('PosStd [m]', fontsize=12, family='Times New Roman')
    ax[0, 0].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax[0, 0].ticklabel_format(style='plain', axis='x')
    ax[0, 0].ticklabel_format(style='plain', axis='y')

    ax[0, 1].plot(time, velstd[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[0, 1].plot(time, velstd[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[0, 1].plot(time, velstd[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[0, 1].grid(True, linestyle='--', color='k', alpha=0.3)
    ax[0, 1].legend(['E', 'N', 'U'], fontsize=12)
    ax[0, 1].set_ylabel('VelStd [m/s]', fontsize=12, family='Times New Roman')
    ax[0, 1].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax[0, 1].ticklabel_format(style='plain', axis='x')
    ax[0, 1].ticklabel_format(style='plain', axis='y')

    ax[1, 0].plot(time, attstd[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[1, 0].plot(time, attstd[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[1, 0].plot(time, attstd[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[1, 0].grid(True, linestyle='--', color='k', alpha=0.3)
    ax[1, 0].legend(['Pitch', 'Roll', 'Yaw'], fontsize=12)
    ax[1, 0].set_ylabel('AttStd [deg]', fontsize=12, family='Times New Roman')
    ax[1, 0].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax[1, 0].ticklabel_format(style='plain', axis='x')
    ax[1, 0].ticklabel_format(style='plain', axis='y')

    ax[2, 0].plot(time, bgstd[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[2, 0].plot(time, bgstd[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[2, 0].plot(time, bgstd[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[2, 0].grid(True, linestyle='--', color='k', alpha=0.3)
    ax[2, 0].legend(['x', 'y', 'z'], fontsize=12)
    ax[2, 0].set_ylabel('BgStd [deg/h]', fontsize=12, family='Times New Roman') 
    ax[2, 0].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax[2, 0].ticklabel_format(style='plain', axis='x')
    ax[2, 0].ticklabel_format(style='plain', axis='y')

    ax[2, 1].plot(time, bastd[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[2, 1].plot(time, bastd[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
    ax[2, 1].plot(time, bastd[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5) 
    ax[2, 1].grid(True, linestyle='--', color='k', alpha=0.3)
    ax[2, 1].legend(['x', 'y', 'z'], fontsize=12)
    ax[2, 1].set_ylabel('BaStd [ug]', fontsize=12, family='Times New Roman')
    ax[2, 1].set_xlabel('GPS Time [s]', fontsize=12, family='Times New Roman')
    ax[2, 1].ticklabel_format(style='plain', axis='x')
    ax[2, 1].ticklabel_format(style='plain', axis='y')

    return fig