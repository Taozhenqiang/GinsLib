import numpy as np
import matplotlib.pyplot as plt
import os
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
                solution[i, 5], solution[i, 6], solution[i, 7], 
                solution[i, 8], solution[i, 9], solution[i, 10]]
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
            # if np.dot(pva_mea[i, 5:8], pva_mea[i, 5:8]) <= 0:
            #     continue
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
            axes[0].set_ylim(-20, 20)

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
            axes[1].set_ylim(-20, 20)

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
            axes[2].set_ylim(-50, 50)

        # Remove scientific notation for axis labels
        for ax in axes:
            ax.ticklabel_format(style='plain', axis='x')
            ax.ticklabel_format(style='plain', axis='y')

        plt.tight_layout()
        figures.append(('pos', fig_pos))

        ## 水平/高程位置误差图 --------------------
        fig_hor_ver, (ax1, ax2) = plt.subplots(2, 1)
        
        # 计算水平位置误差（绝对值）
        horizontal_errors = np.linalg.norm(delta1[:, :2], axis=1)  # E和N方向的2D范数
        
        # 计算高程位置误差（绝对值）
        vertical_errors = np.abs(delta1[:, 2])  # U方向的绝对值
        
        # 绘制水平位置误差
        ax1.plot(t, horizontal_errors, color='#1f77b4', linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        ax1.set_ylabel('Horizontal Error [m]', fontsize=12)
        ax1.grid(True, alpha=0.3)
        ax1.set_title(f'Position Error (GPS week={int(pva_mea[0, 0])}, Epoch={len(pva_mea)})', fontsize=12)
        
        # 添加水平误差统计信息
        hor_max = np.max(horizontal_errors)
        ax1.legend([f'Max: {hor_max:.3f} m'], 
                  loc='upper right', fontsize=10)
        
        # 绘制高程位置误差
        ax2.plot(t, vertical_errors, color='#ff7f0e', linestyle='-', linewidth=1.0, marker='.', markersize=2.5)
        ax2.set_xlabel('GPS Time [s]', fontsize=12)
        ax2.set_ylabel('Vertical Error [m]', fontsize=12)
        ax2.grid(True, alpha=0.3)
        
        # 添加高程误差统计信息
        ver_max = np.max(vertical_errors)
        ax2.legend([f'Max: {ver_max:.3f} m'], 
                  loc='upper right', fontsize=10)
        
        plt.tight_layout()
        figures.append(('3Dpos', fig_hor_ver))

        ## 位置误差CDF图 -----------------
        fig_cdf, (ax1, ax2) = plt.subplots(2,1)
        
        # 计算水平误差的CDF
        sorted_hor_errors = np.sort(horizontal_errors)
        cdf_hor = np.arange(1, len(sorted_hor_errors) + 1) / len(sorted_hor_errors)
        
        # 计算高程误差的CDF
        sorted_ver_errors = np.sort(vertical_errors)
        cdf_ver = np.arange(1, len(sorted_ver_errors) + 1) / len(sorted_ver_errors)
        
        # 绘制水平误差CDF
        ax1.plot(sorted_hor_errors, cdf_hor, color='#1f77b4', linewidth=2)
        
        # 绘制高程误差CDF
        ax2.plot(sorted_ver_errors, cdf_ver, color='#ff7f0e', linewidth=2)
        
        # 标注CEP50/CEP80/CEP95位置
        cep_values = [50, 80, 95]
        cep_colors = ['green', 'orange', 'red']
        
        for cep_val, color in zip(cep_values, cep_colors):
            # 水平误差的CEP值
            hor_cep = calculate_cep(delta1[:, :2], cep_val)
            hor_cdf_idx = np.searchsorted(sorted_hor_errors, hor_cep)
            if hor_cdf_idx < len(cdf_hor):
                hor_cdf_val = cdf_hor[hor_cdf_idx]
                ax1.axvline(x=hor_cep, color=color, linestyle='--', alpha=0.7)
                ax1.plot(hor_cep, hor_cdf_val, 'o', color=color, markersize=8)
                ax1.annotate(f'CEP{cep_val}: {hor_cep:.3f}m', 
                           xy=(hor_cep, hor_cdf_val), 
                           xytext=(10, 10), textcoords='offset points',
                           fontsize=10, color=color,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            
            # 高程误差的CEP值（使用百分位数）
            ver_cep = np.percentile(vertical_errors, cep_val)
            ver_cdf_idx = np.searchsorted(sorted_ver_errors, ver_cep)
            if ver_cdf_idx < len(cdf_ver):
                ver_cdf_val = cdf_ver[ver_cdf_idx]
                ax2.axvline(x=ver_cep, color=color, linestyle='--', alpha=0.7)
                ax2.plot(ver_cep, ver_cdf_val, 's', color=color, markersize=8)
                ax2.annotate(f'CEP{cep_val}: {ver_cep:.3f}m', 
                           xy=(ver_cep, ver_cdf_val), 
                           xytext=(10, -20), textcoords='offset points',
                           fontsize=10, color=color,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        ax1.set_xlabel('Horizontal Error [m]', fontsize=12)
        ax1.set_ylabel('CDF', fontsize=12)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=11)
        # ax1.set_xlim(0, max(np.max(horizontal_errors), np.max(vertical_errors)) * 1.1)
        ax1.set_xlim(0, 10)
        ax1.set_ylim(0, 1.05)
        
        ax2.set_xlabel('Vertical Error [m]', fontsize=12)
        ax2.set_ylabel('CDF', fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=11)
        # ax2.set_xlim(0, max(np.max(horizontal_errors), np.max(vertical_errors)) * 1.1)
        ax2.set_xlim(0, 10)
        ax2.set_ylim(0, 1.05)
        
        plt.tight_layout()
        figures.append(('pos_cdf', fig_cdf))

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

def sync_solutions_by_first(solutions, labels, tolerance=0.01):
    """
    以第一个解决方案为基准，高效同步所有解决方案的时间序列
    使用NumPy向量化操作，避免循环嵌套
    """
    if len(solutions) < 2:
        return solutions, labels
    
    base_solution = solutions[0]
    synced_solutions = [base_solution]
    synced_labels = [labels[0]]
    
    # 向量化计算基准方案的时间序列
    base_times = base_solution[:, 0] * 7 * 24 * 3600 + base_solution[:, 1]
    
    # 同步其他方案
    for i in range(1, len(solutions)):
        solution = solutions[i]
        if solution.shape[0] == 0:
            continue
            
        # 向量化计算当前方案的时间序列
        sol_times = solution[:, 0] * 7 * 24 * 3600 + solution[:, 1]
        
        # 使用广播计算所有时间差矩阵
        # base_times: (m,), sol_times: (n,) -> diff_matrix: (m, n)
        diff_matrix = np.abs(base_times[:, np.newaxis] - sol_times)
        
        # 找到每个基准时间对应的最接近的解决方案时间索引
        min_indices = np.argmin(diff_matrix, axis=1)
        min_diffs = diff_matrix[np.arange(len(base_times)), min_indices]
        
        # 筛选满足容差要求的时间点
        valid_mask = min_diffs < tolerance
        valid_indices = min_indices[valid_mask]
        
        if len(valid_indices) > 0:
            # 提取有效数据
            synced_data = solution[valid_indices]
            
            # 确保时间顺序一致（按基准时间排序）
            sort_order = np.argsort(base_times[valid_mask])
            synced_data = synced_data[sort_order]
            
            synced_solutions.append(synced_data)
            synced_labels.append(labels[i])
    
    return synced_solutions, synced_labels

def plot_err_multi(solutions, labels, reference, flag, time_sync_mode=1):
    """
    多文件误差分析函数：将多个结果文件的误差绘制在同一张图上
    :param solutions: 多个结果文件数据的列表 [solution1, solution2, ...]
    :param labels: 每个结果文件的标签列表 ['方案1', '方案2', ...]
    :param reference: 参考文件数据
    :param flag: 误差类型 ('p', 'v', 'a', 'pv', 'pa', 'va', 'pva')
    :param time_sync_mode: 时间同步模式 (1: 各自独立时间序列; 2: 以第一个文件为基准同步)
    :return: 图形列表, RMS统计字典, CEP统计字典
    """
    
    if len(solutions) != len(labels):
        raise ValueError("解决方案数量和标签数量不匹配")
    
    if len(solutions) == 0:
        raise ValueError("没有提供解决方案数据")

    # 时间同步处理
    if time_sync_mode == 2 and len(solutions) > 1:
        print("使用时间同步模式2：以第一个文件为基准同步时间序列")
        solutions, labels = sync_solutions_by_first(solutions, labels)
    else:
        print("使用时间同步模式1：各自独立时间序列")

    # 初始化图形列表
    figures = []
    rms_stats = {}
    cep_stats = {}
    
    # 处理参考文件
    if reference.shape[0] == 0:
        raise ValueError("参考文件为空")
    
    pva_ref = np.zeros((reference.shape[0], 11))
    for i in range(reference.shape[0]):
        mline = [reference[i, 0], reference[i, 1], reference[i, 2], reference[i, 3], reference[i, 4], 
                reference[i, 5], reference[i, 6], reference[i, 7], 
                reference[i, 8], reference[i, 9], reference[i, 10]]
        pva_ref[i, :] = mline
    
    # 坐标转换
    _, Cne = xyz2blh(pva_ref[0, 2:5])
    
    # 定义颜色列表
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    # ========== 位置误差多文件对比 ==========
    if "p" in flag:
        fig_pos, axes = plt.subplots(3, 1)
        
        # 存储所有解决方案的误差数据
        all_delta_pos = []
        valid_times_list = []  # 新增：存储每个方案的有效时间序列
        
        for i, solution in enumerate(solutions):
            if solution.shape[0] == 0:
                print(f"警告: 解决方案 {labels[i]} 为空，跳过")
                continue
                
            # 将解决方案数据转换为标准格式
            pva_mea = np.zeros((solution.shape[0], 11))
            for j in range(solution.shape[0]):
                mline = [solution[j, 0], solution[j, 1], solution[j, 2], solution[j, 3], solution[j, 4], 
                        solution[j, 5], solution[j, 6], solution[j, 7], 
                        solution[j, 8], solution[j, 9], solution[j, 10]]
                pva_mea[j, :] = mline
            
            # 计算位置误差
            delta_pos = []
            valid_times = []
            
            for j in range(pva_mea.shape[0]):
                if np.dot(pva_mea[j, 2:5], pva_mea[j, 2:5]) <= 0:
                    continue
                    
                time = pva_mea[j, 0] * 7 * 24 * 3600 + pva_mea[j, 1]
                idx = np.where(np.abs(pva_ref[:, 0] * 7 * 24 * 3600 + pva_ref[:, 1] - time) < 0.01)[0]
                
                if len(idx) > 0:
                    if np.dot(pva_ref[idx[0], 2:5], pva_ref[idx[0], 2:5]) <= 0:
                        continue
                    
                    pos_mea = np.dot(Cne, pva_mea[j, 2:5])
                    pos_ref = np.dot(Cne, pva_ref[idx[0], 2:5])
                    error = pos_mea - pos_ref
                    
                    delta_pos.append(error)
                    valid_times.append(pva_mea[j, 1])
            
            if len(delta_pos) > 0:
                delta_pos = np.array(delta_pos)
                all_delta_pos.append(delta_pos)
                valid_times_list.append(valid_times)

                # 绘制东向误差
                axes[0].plot(valid_times, delta_pos[:, 0], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
                
                # 绘制北向误差
                axes[1].plot(valid_times, delta_pos[:, 1], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
                
                # 绘制天向误差
                axes[2].plot(valid_times, delta_pos[:, 2], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])

                # 绘制天向误差
                axes[2].plot(valid_times, delta_pos[:, 2], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
        
        # 计算并存储统计信息
        rms_stats['position_multi'] = {}
        
        # 为每个方案创建包含RMS信息的图例标签
        legend_labels_e = []
        legend_labels_n = []
        legend_labels_u = []

        for i, delta_pos in enumerate(all_delta_pos):
            if len(delta_pos) > 0:
                # RMS计算
                rms_e = np.sqrt(np.sum(delta_pos[:, 0]**2) / delta_pos.shape[0])
                rms_n = np.sqrt(np.sum(delta_pos[:, 1]**2) / delta_pos.shape[0])
                rms_u = np.sqrt(np.sum(delta_pos[:, 2]**2) / delta_pos.shape[0])
                
                # 存储RMS统计信息
                rms_stats['position_multi'][labels[i]] = {
                    'E': rms_e, 'N': rms_n, 'U': rms_u
                }
                
                # 创建包含RMS信息的图例标签
                legend_label_e = f'{labels[i]} (E:{rms_e:.3f}m)'
                legend_labels_e.append(legend_label_e)
                legend_label_n = f'{labels[i]} (N:{rms_n:.3f}m)'
                legend_labels_n.append(legend_label_n)
                legend_label_u = f'{labels[i]} (U:{rms_u:.3f}m)'
                legend_labels_u.append(legend_label_u)

        # 设置子图属性
        axes[0].set_ylabel('E [m]', fontsize=12)
        axes[0].grid(True, alpha=0.3)
        axes[0].legend(fontsize=10)
        axes[0].legend(legend_labels_e, fontsize=9, loc='best')
        
        axes[1].set_ylabel('N [m]', fontsize=12)
        axes[1].grid(True, alpha=0.3)
        axes[1].legend(fontsize=10)
        axes[1].legend(legend_labels_n, fontsize=9, loc='best')
        
        axes[2].set_xlabel('GPS Time [s]', fontsize=12)
        axes[2].set_ylabel('U [m]', fontsize=12)
        axes[2].grid(True, alpha=0.3)
        axes[2].legend(fontsize=10)
        axes[2].legend(legend_labels_u, fontsize=9, loc='best')
        
        plt.tight_layout()
        figures.append(('pos_multi', fig_pos))
        
        # ========== 水平/高程位置误差多文件对比 ==========
        fig_hor_ver, (ax1, ax2) = plt.subplots(2, 1)

        for i, delta_pos in enumerate(all_delta_pos):
            if len(delta_pos) > 0:
                # 计算水平位置误差（绝对值）
                horizontal_errors = np.linalg.norm(delta_pos[:, :2], axis=1)
                # 计算高程位置误差（绝对值）
                vertical_errors = np.abs(delta_pos[:, 2])
                
                # 绘制水平位置误差
                ax1.plot(valid_times_list[i], horizontal_errors, 
                        color=colors[i % len(colors)], linestyle='-', linewidth=1.5, label=labels[i])
                
                # 绘制高程位置误差
                ax2.plot(valid_times_list[i], vertical_errors, 
                        color=colors[i % len(colors)], linestyle='-', linewidth=1.5, label=labels[i])

        # 设置子图属性
        ax1.set_ylabel('Horizontal Error [m]', fontsize=12)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=10)

        ax2.set_xlabel('GPS Time [s]', fontsize=12)
        ax2.set_ylabel('Vertical Error [m]', fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10)

        plt.tight_layout()
        figures.append(('3Dpos_multi', fig_hor_ver))

        # ========== 位置误差CDF多文件对比 ==========
        fig_cdf, (ax1, ax2) = plt.subplots(2, 1)
        
        for i, delta_pos in enumerate(all_delta_pos):
            if len(delta_pos) > 0:
                # 计算水平误差和高程误差
                horizontal_errors = np.linalg.norm(delta_pos[:, :2], axis=1)
                vertical_errors = np.abs(delta_pos[:, 2])
                
                # 水平误差CDF
                sorted_hor_errors = np.sort(horizontal_errors)
                cdf_hor = np.arange(1, len(sorted_hor_errors) + 1) / len(sorted_hor_errors)
                
                # 高程误差CDF
                sorted_ver_errors = np.sort(vertical_errors)
                cdf_ver = np.arange(1, len(sorted_ver_errors) + 1) / len(sorted_ver_errors)
                
                # 绘制水平误差CDF
                ax1.plot(sorted_hor_errors, cdf_hor, color=colors[i % len(colors)], 
                        linewidth=2, label=labels[i])
                
                # 绘制高程误差CDF
                ax2.plot(sorted_ver_errors, cdf_ver, color=colors[i % len(colors)], 
                        linewidth=2, label=labels[i])
                
                # 标注CEP值
                cep_values = [50, 80, 95]
                cep_colors = ['green', 'orange', 'red']
                
                for cep_val, color in zip(cep_values, cep_colors):
                    # 水平误差的CEP值
                    hor_cep = calculate_cep(delta_pos[:, :2], cep_val)
                    hor_cdf_idx = np.searchsorted(sorted_hor_errors, hor_cep)
                    if hor_cdf_idx < len(cdf_hor):
                        hor_cdf_val = cdf_hor[hor_cdf_idx]
                        ax1.axvline(x=hor_cep, color=color, linestyle='--', alpha=0.7)
                        ax1.plot(hor_cep, hor_cdf_val, 'o', color=color, markersize=6)
                        ax1.annotate(f'CEP{cep_val}: {hor_cep:.3f}m', 
                                   xy=(hor_cep, hor_cdf_val), 
                                   xytext=(5, 5), textcoords='offset points',
                                   fontsize=8, color=color,
                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
                    
                    # 高程误差的CEP值（使用百分位数）
                    ver_cep = np.percentile(vertical_errors, cep_val)
                    ver_cdf_idx = np.searchsorted(sorted_ver_errors, ver_cep)
                    if ver_cdf_idx < len(cdf_ver):
                        ver_cdf_val = cdf_ver[ver_cdf_idx]
                        ax2.axvline(x=ver_cep, color=color, linestyle='--', alpha=0.7)
                        ax2.plot(ver_cep, ver_cdf_val, 's', color=color, markersize=6)
                        ax2.annotate(f'CEP{cep_val}: {ver_cep:.3f}m', 
                                   xy=(ver_cep, ver_cdf_val), 
                                   xytext=(5, -15), textcoords='offset points',
                                   fontsize=8, color=color,
                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # 设置子图属性
        ax1.set_xlabel('Horizontal Error [m]', fontsize=12)
        ax1.set_ylabel('CDF', fontsize=12)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=10)
        # ax1.set_xlim(0, max([np.max(np.linalg.norm(delta[:, :2], axis=1)) 
        #                    for delta in all_delta_pos if len(delta) > 0]) * 1.1)
        ax1.set_xlim(0, 20)
        ax1.set_ylim(0, 1.05)
        
        ax2.set_xlabel('Vertical Error [m]', fontsize=12)
        ax2.set_ylabel('CDF', fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10)
        # ax2.set_xlim(0, max([np.max(np.abs(delta[:, 2])) 
        #                    for delta in all_delta_pos if len(delta) > 0]) * 1.1)
        ax2.set_xlim(0, 20)
        ax2.set_ylim(0, 1.05)
        
        plt.tight_layout()
        figures.append(('pos_cdf_multi', fig_cdf))
    
    # ========== 速度误差多文件对比 ==========
    if "v" in flag:
        fig_vel, axes = plt.subplots(3, 1, figsize=(12, 10))
        all_delta_vel = []
        
        for i, solution in enumerate(solutions):
            if solution.shape[0] == 0:
                continue
                
            # 将解决方案数据转换为标准格式
            pva_mea = np.zeros((solution.shape[0], 11))
            for j in range(solution.shape[0]):
                mline = [solution[j, 0], solution[j, 1], solution[j, 2], solution[j, 3], solution[j, 4], 
                        solution[j, 5], solution[j, 6], solution[j, 7], 
                        solution[j, 8], solution[j, 9], solution[j, 10]]
                pva_mea[j, :] = mline
            
            # 计算速度误差
            delta_vel = []
            valid_times = []
            
            for j in range(pva_mea.shape[0]):
                time = pva_mea[j, 0] * 7 * 24 * 3600 + pva_mea[j, 1]
                idx = np.where(np.abs(pva_ref[:, 0] * 7 * 24 * 3600 + pva_ref[:, 1] - time) < 0.01)[0]
                
                if len(idx) > 0:
                    if np.dot(pva_ref[idx[0], 5:8], pva_ref[idx[0], 5:8]) < 0:
                        continue
                    
                    vel_mea = np.dot(Cne, pva_mea[j, 5:8])
                    vel_ref = np.dot(Cne, pva_ref[idx[0], 5:8])
                    error = vel_mea - vel_ref
                    
                    delta_vel.append(error)
                    valid_times.append(pva_mea[j, 1])
            
            if len(delta_vel) > 0:
                delta_vel = np.array(delta_vel)
                all_delta_vel.append(delta_vel)
                
                # 绘制速度误差
                axes[0].plot(valid_times, delta_vel[:, 0], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
                axes[1].plot(valid_times, delta_vel[:, 1], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
                axes[2].plot(valid_times, delta_vel[:, 2], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])

        # 为每个方案创建包含RMS信息的图例标签
        legend_labels_e = []
        legend_labels_n = []
        legend_labels_u = []

        for i, delta_vel in enumerate(all_delta_vel):
            if len(delta_vel) > 0:
                # RMS计算
                rms_e = np.sqrt(np.sum(delta_vel[:, 0]**2) / delta_vel.shape[0])
                rms_n = np.sqrt(np.sum(delta_vel[:, 1]**2) / delta_vel.shape[0])
                rms_u = np.sqrt(np.sum(delta_vel[:, 2]**2) / delta_vel.shape[0])
                
                # 存储RMS统计信息
                rms_stats['position_multi'][labels[i]] = {
                    'E': rms_e, 'N': rms_n, 'U': rms_u
                }
                
                # 创建包含RMS信息的图例标签
                legend_label_e = f'{labels[i]} (E:{rms_e:.3f}m/s)'
                legend_labels_e.append(legend_label_e)
                legend_label_n = f'{labels[i]} (N:{rms_n:.3f}m/s)'
                legend_labels_n.append(legend_label_n)
                legend_label_u = f'{labels[i]} (U:{rms_u:.3f}m/s)'
                legend_labels_u.append(legend_label_u)

        # 设置子图属性
        axes[0].set_ylabel('E [m/s]', fontsize=12)
        axes[0].grid(True, alpha=0.3)
        axes[0].legend(fontsize=10)
        axes[0].legend(legend_labels_e, fontsize=9, loc='best')
        
        axes[1].set_ylabel('N [m/s]', fontsize=12)
        axes[1].grid(True, alpha=0.3)
        axes[1].legend(fontsize=10)
        axes[1].legend(legend_labels_n, fontsize=9, loc='best')
        
        axes[2].set_xlabel('GPS Time [s]', fontsize=12)
        axes[2].set_ylabel('U [m/s]', fontsize=12)
        axes[2].grid(True, alpha=0.3)
        axes[2].legend(fontsize=10)
        axes[2].legend(legend_labels_u, fontsize=9, loc='best')
        
        plt.tight_layout()
        figures.append(('vel_multi', fig_vel))
    
    # ========== 姿态误差多文件对比 ==========
    if "a" in flag:
        fig_att, axes = plt.subplots(3, 1, figsize=(12, 10))
        all_delta_att = []
        
        for i, solution in enumerate(solutions):
            if solution.shape[0] == 0:
                continue
                
            # 将解决方案数据转换为标准格式
            pva_mea = np.zeros((solution.shape[0], 11))
            for j in range(solution.shape[0]):
                mline = [solution[j, 0], solution[j, 1], solution[j, 2], solution[j, 3], solution[j, 4], 
                        solution[j, 5], solution[j, 6], solution[j, 7], 
                        solution[j, 8], solution[j, 9], solution[j, 10]]
                pva_mea[j, :] = mline
            
            # 计算姿态误差
            delta_att = []
            valid_times = []
            
            for j in range(pva_mea.shape[0]):
                if np.dot(pva_mea[j, 8:11], pva_mea[j, 8:11]) <= 0:
                    continue
                    
                time = pva_mea[j, 0] * 7 * 24 * 3600 + pva_mea[j, 1]
                idx = np.where(np.abs(pva_ref[:, 0] * 7 * 24 * 3600 + pva_ref[:, 1] - time) < 0.01)[0]
                
                if len(idx) > 0:
                    if np.dot(pva_ref[idx[0], 8:11], pva_ref[idx[0], 8:11]) <= 0:
                        continue
                    
                    att_mea = pva_mea[j, 8:11]
                    att_ref = pva_ref[idx[0], 8:11]
                    error = att_mea - att_ref
                    
                    # 处理航向角跳变
                    if error[2] >= 300:
                        error[2] -= 360
                    elif error[2] <= -300:
                        error[2] += 360
                    
                    delta_att.append(error)
                    valid_times.append(pva_mea[j, 1])
            
            if len(delta_att) > 0:
                delta_att = np.array(delta_att)
                
                # 绘制姿态误差
                axes[0].plot(valid_times, delta_att[:, 0], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
                axes[1].plot(valid_times, delta_att[:, 1], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
                axes[2].plot(valid_times, delta_att[:, 2], color=colors[i % len(colors)], 
                           linestyle='-', linewidth=1.0, marker='.', markersize=2.5, label=labels[i])
        

        # 为每个方案创建包含RMS信息的图例标签
        legend_labels_pitch = []
        legend_labels_roll = []
        legend_labels_heading = []

        for i, delta_att in enumerate(all_delta_att):
            if len(delta_att) > 0:
                # RMS计算
                rms_e = np.sqrt(np.sum(delta_att[:, 0]**2) / delta_att.shape[0])
                rms_n = np.sqrt(np.sum(delta_att[:, 1]**2) / delta_att.shape[0])
                rms_u = np.sqrt(np.sum(delta_att[:, 2]**2) / delta_att.shape[0])
                
                # 存储RMS统计信息
                rms_stats['position_multi'][labels[i]] = {
                    'E': rms_e, 'N': rms_n, 'U': rms_u
                }
                
                # 创建包含RMS信息的图例标签
                legend_label_pitch = f'{labels[i]} (E:{rms_e:.3f}°)'
                legend_labels_pitch.append(legend_label_pitch)
                legend_label_roll = f'{labels[i]} (N:{rms_n:.3f}°)'   
                legend_labels_roll.append(legend_label_roll)
                legend_label_heading = f'{labels[i]} (U:{rms_u:.3f}°)'
                legend_labels_heading.append(legend_label_heading)

        # 设置子图属性
        axes[0].set_ylabel('Pitch [°]', fontsize=12)
        axes[0].grid(True, alpha=0.3)
        axes[0].legend(fontsize=10)
        axes[0].legend(legend_labels_pitch, fontsize=9, loc='best')
        
        axes[1].set_ylabel('Roll [°]', fontsize=12)
        axes[1].grid(True, alpha=0.3)
        axes[1].legend(fontsize=10)
        axes[1].legend(legend_labels_roll, fontsize=9, loc='best')
        
        axes[2].set_xlabel('GPS Time [s]', fontsize=12)
        axes[2].set_ylabel('Heading [°]', fontsize=12)
        axes[2].grid(True, alpha=0.3)
        axes[2].legend(fontsize=10)
        axes[2].legend(legend_labels_heading, fontsize=9, loc='best')
        
        plt.tight_layout()
        figures.append(('att_multi', fig_att))
    
    return figures, rms_stats, cep_stats