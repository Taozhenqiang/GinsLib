import numpy as np
import matplotlib.pyplot as plt
from xyz2blh import xyz2blh

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

    plt.plot(pos0[:, 0] - pos0[0, 0], pos0[:, 1] - pos0[0, 1], ':', linewidth=0.1, color=[0.5, 0.5, 0.5])
    plt.plot(pos0[:, 0] - pos0[0, 0], pos0[:, 1] - pos0[0, 1], '.b')

    plt.grid(True, linestyle='--', color='k', alpha=0.3)
    plt.xlabel('E [m]', fontsize=12, family='Times New Roman')
    plt.ylabel('N [m]', fontsize=12, family='Times New Roman')
    plt.title('Trajectory', fontsize=12, family='Times New Roman')
    plt.axis('equal')
    plt.gca().ticklabel_format(style='plain', axis='x')  # Disable scientific notation for x-axis

    plt.show()

