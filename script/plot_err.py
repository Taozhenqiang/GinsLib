import numpy as np
import matplotlib.pyplot as plt
from xyz2blh import xyz2blh


def plot_err(solution, reference, flag):

    nsol = solution.shape[0]
    nref = reference.shape[0]
    
    if nsol == 0:
        raise ValueError("Solution is empty!!!")
    if nref == 0:
        raise ValueError("Reference truth is empty!!!")
    
    pva_mea = np.zeros((nsol, 11))
    pva_ref = np.zeros((nref, 11))

    for i in range(nsol):
        mline= [solution[i, 0], solution[i, 1], solution[i, 2], solution[i, 3], solution[i, 4], 
                solution[i, 6], solution[i, 7], solution[i, 8], 
                solution[i, 9], solution[i, 10], solution[i, 11]]
        pva_mea[i, :] = mline

    for i in range(nref):
        mline = [reference[i, 0], reference[i, 1], reference[i, 2], reference[i, 3], reference[i, 4], 
                reference[i, 5], reference[i, 6], reference[i, 7], 
                reference[i, 8], reference[i, 9], reference[i, 10]]
        pva_ref[i, :] = mline

    # 坐标转换
    _, Cne = xyz2blh(pva_ref[0, 2:5])
    
    t = np.zeros(nsol)
    pos1 = np.zeros((nsol, 3))
    pos2 = np.zeros((nsol, 3))
    vel1 = np.zeros((nsol, 3))
    vel2 = np.zeros((nsol, 3))
    datt = np.zeros((nsol, 3))

    if "p" in flag:
        k = 0
        for i in range(nsol):
            if np.dot(pva_mea[i, 2:5], pva_mea[i, 2:5]) <= 0:
                continue
            time = pva_mea[i, 0] * 7 * 24 * 3600 + pva_mea[i, 1]
            idx = np.where(np.abs(pva_ref[:, 0] * 7 * 24 * 3600 + pva_ref[:, 1] - time) < 0.01)[0]
            if len(idx) > 0:
                if np.dot(pva_ref[idx[0], 2:5], pva_ref[idx[0], 2:5]) <= 0:
                    continue
                t[k] = pva_mea[i, 1]
                pos1[k, :] = np.dot(Cne, pva_mea[i, 2:5])
                pos2[k, :] = np.dot(Cne, pva_ref[idx[0], 2:5])
                k += 1
        if k == 0:
            raise ValueError("The position solution or reference truth does not exist!!!")
        if k < nsol:
            t = t[:k]
            pos1 = pos1[:k, :]
            pos2 = pos2[:k, :]

    if "v" in flag:
        k = 0
        for i in range(nsol):
            if np.dot(pva_mea[i, 5:8], pva_mea[i, 5:8]) <= 0:
                continue
            time = pva_mea[i, 0] * 7 * 24 * 3600 + pva_mea[i, 1]
            idx = np.where(np.abs(pva_ref[:, 0] * 7 * 24 * 3600 + pva_ref[:, 1] - time) < 0.01)[0]
            if len(idx) > 0:
                if np.dot(pva_ref[idx[0], 5:8], pva_ref[idx[0], 5:8]) < 0:
                    continue
                t[k] = pva_mea[i, 1]
                vel1[k, :] = np.dot(Cne, pva_mea[i, 5:8])
                vel2[k, :] = np.dot(Cne, pva_ref[idx[0], 5:8])
                k += 1
        if k == 0:
            raise ValueError("The velocity solution or reference truth does not exist!!!")
        if k < nsol:
            t = t[:k]
            vel1 = vel1[:k, :]
            vel2 = vel2[:k, :]

    if "a" in flag:
        k = 0
        for i in range(nsol):
            if np.dot(pva_mea[i, 8:11], pva_mea[i, 8:11]) <= 0:
                continue
            time = pva_mea[i, 0] * 7 * 24 * 3600 + pva_mea[i, 1]
            idx = np.where(np.abs(pva_ref[:, 0] * 7 * 24 * 3600 + pva_ref[:, 1] - time) < 0.01)[0]
            if len(idx) > 0:
                if np.dot(pva_ref[idx[0], 8:11], pva_ref[idx[0], 8:11]) <= 0:
                    continue
                t[k] = pva_mea[i, 1]
                att1 = pva_mea[i, 8:11]
                att2 = pva_ref[idx[0], 8:11]
                datt[k, :] = att1 - att2
                k += 1
        if k == 0:
            raise ValueError("The attitude solution or reference truth does not exist!!!")
        if k < nsol:
            t = t[:k]
            datt = datt[:k, :]

    # 位置误差
    if "p" in flag:
        delta1 = pos1 - pos2

        fig, axes = plt.subplots(3, 1)
        Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

        # 剔除异常值
        # idx = np.linalg.norm(delta1, axis=1) / np.linalg.norm(np.mean(np.abs(delta1), axis=0)) > 100
        # delta11 = delta1[~idx, :]
        delta11 = delta1
        
        # 绘制位置误差
        axes[0].plot(t, delta1[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[0].set_ylabel('E [m]')
        axes[0].grid(True)
        axes[0].set_title(f"Position error (GPS week={int(pva_mea[0, 0])})")
        axes[0].legend([f'RMS: {np.sqrt(np.sum(delta11[:, 0]**2) / delta11.shape[0]):.4f} m'])

        axes[1].plot(t, delta1[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[1].set_ylabel('N [m]')
        axes[1].grid(True)
        axes[1].legend([f'RMS: {np.sqrt(np.sum(delta11[:, 1]**2) / delta11.shape[0]):.4f} m'])

        axes[2].plot(t, delta1[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[2].set_xlabel('GPS Time [s]')
        axes[2].set_ylabel('U [m]')
        axes[2].grid(True)
        axes[2].legend([f'RMS: {np.sqrt(np.sum(delta11[:, 2]**2) / delta11.shape[0]):.4f} m'])

        # Remove scientific notation for axis labels
        for ax in axes:
            ax.ticklabel_format(style='plain', axis='x')
            ax.ticklabel_format(style='plain', axis='y')

    if "v" in flag:
        delta2 = vel1 - vel2

        fig, axes = plt.subplots(3, 1)
        Fcolor = ["#ffcc66", "#14a959", "#ff6666"]
        
        # 绘制位置误差
        axes[0].plot(t, delta2[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[0].set_ylabel('E [m/s]')
        axes[0].grid(True)
        axes[0].set_title(f"Velocity error (GPS week={int(pva_mea[0, 0])})")
        axes[0].legend([f'RMS: {np.sqrt(np.sum(delta2[:, 0]**2) / delta2.shape[0]):.4f} m/s'])

        axes[1].plot(t, delta2[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[1].set_ylabel('N [m/s]')
        axes[1].grid(True)
        axes[1].legend([f'RMS: {np.sqrt(np.sum(delta2[:, 1]**2) / delta2.shape[0]):.4f} m/s'])

        axes[2].plot(t, delta2[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[2].set_xlabel('GPS Time [s]')
        axes[2].set_ylabel('U [m/s]')
        axes[2].grid(True)
        axes[2].legend([f'RMS: {np.sqrt(np.sum(delta2[:, 2]**2) / delta2.shape[0]):.4f} m/s'])

        # Remove scientific notation for axis labels
        for ax in axes:
            ax.ticklabel_format(style='plain', axis='x')
            ax.ticklabel_format(style='plain', axis='y')

    if "a" in flag:
        delta3 = datt

        # 条件1：第三列大于或等于 300
        idx_att = delta3[:, 2] >= 300
        delta3[idx_att, 2] -= 360  # 减去 360

        # 条件2：第三列小于或等于 -300
        idx_att = delta3[:, 2] <= -300
        delta3[idx_att, 2] += 360  # 加上 360

        fig, axes = plt.subplots(3, 1)
        Fcolor = ["#ffcc66", "#14a959", "#ff6666"]
        
        # 绘制位置误差
        axes[0].plot(t, delta3[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[0].set_ylabel('Pitch [deg]')
        axes[0].grid(True)
        axes[0].set_title(f"Attitude error (GPS week={int(pva_mea[0, 0])})")
        axes[0].legend([f'RMS: {np.sqrt(np.sum(delta3[:, 0]**2) / delta3.shape[0]):.4f} deg'])

        axes[1].plot(t, delta3[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[1].set_ylabel('Roll [deg]')
        axes[1].grid(True)
        axes[1].legend([f'RMS: {np.sqrt(np.sum(delta3[:, 1]**2) / delta3.shape[0]):.4f} deg'])

        axes[2].plot(t, delta3[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[2].set_xlabel('GPS Time [s]')
        axes[2].set_ylabel('Yaw[deg]')
        axes[2].grid(True)
        axes[2].legend([f'RMS: {np.sqrt(np.sum(delta3[:, 2]**2) / delta3.shape[0]):.4f} deg'])

        # Remove scientific notation for axis labels
        for ax in axes:
            ax.ticklabel_format(style='plain', axis='x')
            ax.ticklabel_format(style='plain', axis='y')            


    # Adjust spacing
    plt.tight_layout()

    # Show the plot
    plt.show()