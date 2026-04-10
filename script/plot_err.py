import numpy as np
import matplotlib.pyplot as plt
from coordinate_transformation import xyz2blh

def calculate_cep(errors, percentage):
    """
    计算给定百分比的CEP (Circular Error Probable)
    :param errors: 二维误差数组 (N x 2)
    :param percentage: 百分比 (0-100)
    :return: CEP值
    """
    # 计算每个点的水平误差
    horizontal_errors = np.linalg.norm(errors, axis=1)
    # 按升序排序
    sorted_errors = np.sort(horizontal_errors)
    # 计算对应百分比的索引
    idx = int(len(sorted_errors) * percentage / 100)
    # 返回对应索引的值
    return sorted_errors[idx] if idx < len(sorted_errors) else sorted_errors[-1]

def plot_err(solution, reference, flag, manual_yaxis=False):

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

    # 初始化RMS统计字典
    rms_stats = {}
    # 初始化CEP统计字典
    cep_stats = {}
    
    # 初始化图形列表
    figures = []

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

        fig_pos, axes = plt.subplots(3, 1)
        Fcolor = ["#ffcc66", "#14a959", "#ff6666"]

        # 剔除异常值
        # idx = np.linalg.norm(delta1, axis=1) / np.linalg.norm(np.mean(np.abs(delta1), axis=0)) > 100
        # delta11 = delta1[~idx, :]
        delta11 = delta1
        
        # 计算位置误差RMS值
        rms_e = np.sqrt(np.sum(delta11[:, 0]**2) / delta11.shape[0])
        rms_n = np.sqrt(np.sum(delta11[:, 1]**2) / delta11.shape[0])
        rms_u = np.sqrt(np.sum(delta11[:, 2]**2) / delta11.shape[0])
        rms_3d = np.sqrt(np.sum(np.linalg.norm(delta11, axis=1)**2) / delta11.shape[0])      
        
        # 存储位置误差RMS值
        rms_stats['position'] = {
            'E': rms_e,
            'N': rms_n, 
            'U': rms_u,
            '3D': rms_3d
        }

        # 计算位置误差CEP值
        # 水平方向（东和北）
        horizontal_errors = delta11[:, :2]  # 只取东和北方向
        cep50_horizontal = calculate_cep(horizontal_errors, 50)
        cep68_horizontal = calculate_cep(horizontal_errors, 68)
        cep80_horizontal = calculate_cep(horizontal_errors, 80)
        cep95_horizontal = calculate_cep(horizontal_errors, 95)
        cep99_horizontal = calculate_cep(horizontal_errors, 99)
        
        # 高程方向（天向）
        vertical_errors = delta11[:, 2:3]  # 只取天向
        cep50_vertical = np.percentile(np.abs(vertical_errors), 50)
        cep68_vertical = np.percentile(np.abs(vertical_errors), 68)
        cep80_vertical = np.percentile(np.abs(vertical_errors), 80)
        cep95_vertical = np.percentile(np.abs(vertical_errors), 95)
        cep99_vertical = np.percentile(np.abs(vertical_errors), 99)
        
        # 存储位置误差CEP值
        cep_stats['position'] = {
            'horizontal': {
                'CEP50': cep50_horizontal,
                'CEP68': cep68_horizontal,
                'CEP80': cep80_horizontal,
                'CEP95': cep95_horizontal,
                'CEP99': cep99_horizontal
            },
            'vertical': {
                'CEP50': cep50_vertical,
                'CEP68': cep68_vertical,
                'CEP80': cep80_vertical,
                'CEP95': cep95_vertical,
                'CEP99': cep99_vertical
            }
        }
        
        # 绘制位置误差
        axes[0].plot(t, delta1[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[0].set_ylabel('E [m]')
        axes[0].grid(True)
        axes[0].set_title(f"Position error (GPS week={int(pva_mea[0, 0])})")
        axes[0].legend([f'RMS: {rms_e:.4f} m'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            # 计算y轴范围，确保包含所有数据点
            y_min = np.min(delta1[:, 0])
            y_max = np.max(delta1[:, 0])
            y_range = y_max - y_min
            # 设置y轴范围，留出10%的边距
            # axes[0].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
            axes[0].set_ylim(-50, 50)

        axes[1].plot(t, delta1[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[1].set_ylabel('N [m]')
        axes[1].grid(True)
        axes[1].legend([f'RMS: {rms_n:.4f} m'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta1[:, 1])
            y_max = np.max(delta1[:, 1])
            y_range = y_max - y_min
            # axes[1].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
            axes[1].set_ylim(-50, 50)

        axes[2].plot(t, delta1[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[2].set_xlabel('GPS Time [s]')
        axes[2].set_ylabel('U [m]')
        axes[2].grid(True)
        axes[2].legend([f'RMS: {rms_u:.4f} m'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta1[:, 2])
            y_max = np.max(delta1[:, 2])
            y_range = y_max - y_min
            # axes[2].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
            axes[2].set_ylim(-100, 100)

        # Remove scientific notation for axis labels
        for ax in axes:
            ax.ticklabel_format(style='plain', axis='x')
            ax.ticklabel_format(style='plain', axis='y')

        plt.tight_layout()
        figures.append(('pos', fig_pos))

    if "v" in flag:
        delta2 = vel1 - vel2

        fig_vel, axes = plt.subplots(3, 1)
        Fcolor = ["#ffcc66", "#14a959", "#ff6666"]
        
        # 计算速度误差RMS值
        rms_e = np.sqrt(np.sum(delta2[:, 0]**2) / delta2.shape[0])
        rms_n = np.sqrt(np.sum(delta2[:, 1]**2) / delta2.shape[0])
        rms_u = np.sqrt(np.sum(delta2[:, 2]**2) / delta2.shape[0])
        rms_3d = np.sqrt(np.sum(np.linalg.norm(delta2, axis=1)**2) / delta2.shape[0])
        
        # 存储速度误差RMS值
        rms_stats['velocity'] = {
            'E': rms_e,
            'N': rms_n,
            'U': rms_u,
            '3D': rms_3d
        }
        # 计算速度误差CEP值
        # 水平方向（东和北）
        horizontal_errors = delta2[:, :2]  # 只取东和北方向
        cep50_horizontal = calculate_cep(horizontal_errors, 50)
        cep68_horizontal = calculate_cep(horizontal_errors, 68)
        cep80_horizontal = calculate_cep(horizontal_errors, 80)
        cep95_horizontal = calculate_cep(horizontal_errors, 95)
        cep99_horizontal = calculate_cep(horizontal_errors, 99)
        
        # 高程方向（天向）
        vertical_errors = delta2[:, 2:3]  # 只取天向
        cep50_vertical = np.percentile(np.abs(vertical_errors), 50)
        cep68_vertical = np.percentile(np.abs(vertical_errors), 68)
        cep80_vertical = np.percentile(np.abs(vertical_errors), 80)
        cep95_vertical = np.percentile(np.abs(vertical_errors), 95)
        cep99_vertical = np.percentile(np.abs(vertical_errors), 99)
        
        # 存储速度误差CEP值
        cep_stats['velocity'] = {
            'horizontal': {
                'CEP50': cep50_horizontal,
                'CEP68': cep68_horizontal,
                'CEP80': cep80_horizontal,
                'CEP95': cep95_horizontal,
                'CEP99': cep99_horizontal
            },
            'vertical': {
                'CEP50': cep50_vertical,
                'CEP68': cep68_vertical,
                'CEP80': cep80_vertical,
                'CEP95': cep95_vertical,
                'CEP99': cep99_vertical
            }
        }
        
        # 绘制速度误差
        axes[0].plot(t, delta2[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[0].set_ylabel('E [m/s]')
        axes[0].grid(True)
        axes[0].set_title(f"Velocity error (GPS week={int(pva_mea[0, 0])})")
        axes[0].legend([f'RMS: {rms_e:.4f} m/s'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta2[:, 0])
            y_max = np.max(delta2[:, 0])
            y_range = y_max - y_min
            axes[0].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        axes[1].plot(t, delta2[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[1].set_ylabel('N [m/s]')
        axes[1].grid(True)
        axes[1].legend([f'RMS: {rms_n:.4f} m/s'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta2[:, 1])
            y_max = np.max(delta2[:, 1])
            y_range = y_max - y_min
            axes[1].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        axes[2].plot(t, delta2[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[2].set_xlabel('GPS Time [s]')
        axes[2].set_ylabel('U [m/s]')
        axes[2].grid(True)
        axes[2].legend([f'RMS: {rms_u:.4f} m/s'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta2[:, 2])
            y_max = np.max(delta2[:, 2])
            y_range = y_max - y_min
            axes[2].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        # Remove scientific notation for axis labels
        for ax in axes:
            ax.ticklabel_format(style='plain', axis='x')
            ax.ticklabel_format(style='plain', axis='y')

        plt.tight_layout()
        figures.append(('vel', fig_vel))

    if "a" in flag:
        delta3 = datt

        # 条件1：第三列大于或等于 300
        idx_att = delta3[:, 2] >= 300
        delta3[idx_att, 2] -= 360  # 减去 360

        # 条件2：第三列小于或等于 -300
        idx_att = delta3[:, 2] <= -300
        delta3[idx_att, 2] += 360  # 加上 360

        fig_att, axes = plt.subplots(3, 1)
        Fcolor = ["#ffcc66", "#14a959", "#ff6666"]
        
        # 计算姿态误差RMS值
        rms_pitch = np.sqrt(np.sum(delta3[:, 0]**2) / delta3.shape[0])
        rms_roll = np.sqrt(np.sum(delta3[:, 1]**2) / delta3.shape[0])
        rms_yaw = np.sqrt(np.sum(delta3[:, 2]**2) / delta3.shape[0])
        
        # 存储姿态误差RMS值
        rms_stats['attitude'] = {
            'pitch': rms_pitch,
            'roll': rms_roll,
            'yaw': rms_yaw
        }
        # 计算姿态误差CEP值
        # 对于姿态，我们计算每个轴的CEP50、CEP68和CEP95
        pitch_cep50 = np.percentile(np.abs(delta3[:, 0]), 50)
        pitch_cep68 = np.percentile(np.abs(delta3[:, 0]), 68)
        pitch_cep80 = np.percentile(np.abs(delta3[:, 0]), 80)
        pitch_cep95 = np.percentile(np.abs(delta3[:, 0]), 95)
        pitch_cep99 = np.percentile(np.abs(delta3[:, 0]), 99)
        
        roll_cep50 = np.percentile(np.abs(delta3[:, 1]), 50)
        roll_cep68 = np.percentile(np.abs(delta3[:, 1]), 68)
        roll_cep80 = np.percentile(np.abs(delta3[:, 1]), 80)
        roll_cep95 = np.percentile(np.abs(delta3[:, 1]), 95)
        roll_cep99 = np.percentile(np.abs(delta3[:, 1]), 99)
        
        yaw_cep50 = np.percentile(np.abs(delta3[:, 2]), 50)
        yaw_cep68 = np.percentile(np.abs(delta3[:, 2]), 68)
        yaw_cep80 = np.percentile(np.abs(delta3[:, 2]), 80)
        yaw_cep95 = np.percentile(np.abs(delta3[:, 2]), 95)
        yaw_cep99 = np.percentile(np.abs(delta3[:, 2]), 99)
        
        # 存储姿态误差CEP值
        cep_stats['attitude'] = {
            'pitch': {
                'CEP50': pitch_cep50,
                'CEP68': pitch_cep68,
                'CEP80': pitch_cep80,
                'CEP95': pitch_cep95,
                'CEP99': pitch_cep99
            },
            'roll': {
                'CEP50': roll_cep50,
                'CEP68': roll_cep68,
                'CEP80': roll_cep80,
                'CEP95': roll_cep95,
                'CEP99': roll_cep99
            },
            'yaw': {
                'CEP50': yaw_cep50,
                'CEP68': yaw_cep68,
                'CEP80': yaw_cep80,
                'CEP95': yaw_cep95,
                'CEP99': yaw_cep99
            }
        }
        
        # 绘制姿态误差
        axes[0].plot(t, delta3[:, 0], color=Fcolor[0], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[0].set_ylabel('Pitch [deg]')
        axes[0].grid(True)
        axes[0].set_title(f"Attitude error (GPS week={int(pva_mea[0, 0])})")
        axes[0].legend([f'RMS: {rms_pitch:.4f} deg'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta3[:, 0])
            y_max = np.max(delta3[:, 0])
            y_range = y_max - y_min
            axes[0].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        axes[1].plot(t, delta3[:, 1], color=Fcolor[1], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[1].set_ylabel('Roll [deg]')
        axes[1].grid(True)
        axes[1].legend([f'RMS: {rms_roll:.4f} deg'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta3[:, 1])
            y_max = np.max(delta3[:, 1])
            y_range = y_max - y_min
            axes[1].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        axes[2].plot(t, delta3[:, 2], color=Fcolor[2], linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        axes[2].set_xlabel('GPS Time [s]')
        axes[2].set_ylabel('Yaw[deg]')
        axes[2].grid(True)
        axes[2].legend([f'RMS: {rms_yaw:.4f} deg'])
        
        # 如果启用手动设置y轴标签
        if manual_yaxis:
            y_min = np.min(delta3[:, 2])
            y_max = np.max(delta3[:, 2])
            y_range = y_max - y_min
            axes[2].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        # Remove scientific notation for axis labels
        for ax in axes:
            ax.ticklabel_format(style='plain', axis='x')
            ax.ticklabel_format(style='plain', axis='y')

        plt.tight_layout()
        figures.append(('att', fig_att))

    return figures, rms_stats, cep_stats