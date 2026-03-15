import os
import numpy as np
import logging
from read_sol import read_solution
from read_ref import read_ref
from plot_solution import plot_trajectory
from plot_solution import plot_position
from plot_solution import plot_velocity
from plot_solution import plot_attitude
from plot_solution import plot_imubias
from plot_err import plot_err
import matplotlib.pyplot as plt

def batch_plot_analysis(sol_file_path, ref_file_path, process_mode):
    """
    批处理解算可视化分析函数
    
    参数:
    - sol_file_path: 结果文件路径
    - ref_file_path: 参考文件路径  
    - process_mode: 解算模式字符串（包含LC/TC表示GNSS/INS，否则为GNSS）
    
    返回:
    - error_stats: 误差统计信息字典
    """
    
    # 提取data_name：根据处理模式提取不同的部分
    file_name = os.path.basename(sol_file_path)
    parts = file_name.split('_')
    
    if 'LC' in process_mode or 'TC' in process_mode:
        # GNSS/INS模式：提取倒数第三个_前面的部分
        if len(parts) >= 3:
            data_name = '_'.join(parts[:-2])  # 取除了最后两个部分的所有部分
        else:
            data_name = file_name
    else:
        # GNSS模式：提取倒数第二个_前面的部分
        if len(parts) >= 2:
            data_name = '_'.join(parts[:-1])  # 取除了最后一个部分的所有部分
        else:
            data_name = file_name

    # 提取data_name2：sol_file_path中最后一个分隔符后面的部分
    data_name2 = os.path.basename(sol_file_path)

    # 去掉.pos后缀，用于文件名构建
    if data_name2.endswith('.pos'):
        data_name_without_ext = data_name2[:-4]  # 去掉最后4个字符(.pos)
    else:
        data_name_without_ext = data_name2
    
    # 初始化误差统计信息
    error_stats = {
        'data_name': data_name,
        'process_mode': process_mode,
        'position_rms': None,
        'velocity_rms': None,
        'attitude_rms': None,
        'position_cep': None,
        'velocity_cep': None,
        'attitude_cep': None
    }
    
    logging.info(f"{'='*60}")
    logging.info(f"开始处理: {error_stats['data_name']}")
    logging.info(f"处理模式: {process_mode}")
    logging.info(f"{'='*60}")
    
    # 根据处理模式确定绘制选项
    if 'LC' in process_mode or 'TC' in process_mode:
        # GNSS/INS模式：绘制IMU零偏、位置、速度与姿态误差
        plot_options = ['bias', 'err']
        error_type = 'pva'  # 位置、速度、姿态误差
        logging.info("GNSS/INS模式，绘制IMU零偏和位置、速度、姿态误差")
    else:
        # GNSS模式：仅绘制位置误差
        plot_options = ['err']
        error_type = 'p'  # 仅位置误差
        logging.info("GNSS模式，仅绘制位置误差")
    
    # 结果文件读取配置
    skip_lines = 14
    row = 10000
    # GINLIB: gps week, sow, pos[x/y/z], ratio, vel[x/y/z], att[pitch/roll/heading], bg[x/y/z], ba[x/y/z]
    col = 18
    sample = 1/100

    # 读取结果文件
    data = read_solution(sol_file_path, skip_lines, row, col, sample)

    # 检查是否成功读取
    if data is not None:
        logging.info(f"✓ 结果文件数据读取成功 {sol_file_path}")
    else:
        logging.error(f"✗ 结果文件数据读取失败 {sol_file_path}")
        return error_stats

    # 参考文件读取配置
    skip_lines = 0
    row = 10000
    # gps week, sow, pos[x/y/z], vel[x/y/z], att[pitch/roll/heading]
    col = 11
    ins_flag = 1
    ins_interval = 1/100
    GNSS_interval = 1
    # truth.truth
    idx = [0,1,2,3,4,5,6,7,8,9,10]

    # 读取参考文件
    ref_data = read_ref(ref_file_path, skip_lines, row, col, ins_flag, ins_interval, GNSS_interval, idx)

    # 检查是否成功读取
    if ref_data is not None:
        logging.info(f"✓ 参考文件数据读取成功 {ref_file_path}")
    else:
        logging.error(f"✗ 参考文件数据读取失败 {ref_file_path}")   
        return error_stats
    
    figures = []  # 存储所有图形对象和对应的类型

    for plot_type in plot_options:
        if plot_type == 'bias':
            logging.info("绘制IMU零偏图...")
            # plot_imubias现在返回两个图形对象
            fig_bg, fig_ba = plot_imubias(data)
            figures.append(('bg', fig_bg))  # 陀螺仪零偏图
            figures.append(('ba', fig_ba))  # 加速度计零偏图
        elif plot_type == 'err':
            # plot_err现在返回图形列表和RMS统计信息
            err_figures, rms_stats, cep_stats = plot_err(data, ref_data, error_type, False)
            figures.extend(err_figures)  # 添加所有误差图
            
            # 将rms_stats信息赋值给error_stats
            if 'position' in rms_stats:
                error_stats['position_rms'] = rms_stats['position']
                logging.info(f"位置误差RMS: E={rms_stats['position']['E']:.4f}m, N={rms_stats['position']['N']:.4f}m, U={rms_stats['position']['U']:.4f}m, 3D={rms_stats['position']['3D']:.4f}m")
            
            if 'velocity' in rms_stats:
                error_stats['velocity_rms'] = rms_stats['velocity']
                logging.info(f"速度误差RMS: E={rms_stats['velocity']['E']:.4f}m/s, N={rms_stats['velocity']['N']:.4f}m/s, U={rms_stats['velocity']['U']:.4f}m/s, 3D={rms_stats['velocity']['3D']:.4f}m/s")
            
            if 'attitude' in rms_stats:
                error_stats['attitude_rms'] = rms_stats['attitude']
                logging.info(f"姿态误差RMS: Pitch={rms_stats['attitude']['pitch']:.4f}°, Roll={rms_stats['attitude']['roll']:.4f}°, Yaw={rms_stats['attitude']['yaw']:.4f}°")
            
            # 将cep_stats信息赋值给error_stats
            if 'position' in cep_stats:
                error_stats['position_cep'] = cep_stats['position']
                logging.info(f"位置误差CEP: 水平CEP50={cep_stats['position']['horizontal']['CEP50']:.4f}m, 水平CEP68={cep_stats['position']['horizontal']['CEP68']:.4f}m, 水平CEP95={cep_stats['position']['horizontal']['CEP95']:.4f}m")
                logging.info(f"位置误差CEP: 垂直CEP50={cep_stats['position']['vertical']['CEP50']:.4f}m, 垂直CEP68={cep_stats['position']['vertical']['CEP68']:.4f}m, 垂直CEP95={cep_stats['position']['vertical']['CEP95']:.4f}m")
            
            if 'velocity' in cep_stats:
                error_stats['velocity_cep'] = cep_stats['velocity']
                logging.info(f"速度误差CEP: 水平CEP50={cep_stats['velocity']['horizontal']['CEP50']:.4f}m/s, 水平CEP68={cep_stats['velocity']['horizontal']['CEP68']:.4f}m/s, 水平CEP95={cep_stats['velocity']['horizontal']['CEP95']:.4f}m/s")
                logging.info(f"速度误差CEP: 垂直CEP50={cep_stats['velocity']['vertical']['CEP50']:.4f}m/s, 垂直CEP68={cep_stats['velocity']['vertical']['CEP68']:.4f}m/s, 垂直CEP95={cep_stats['velocity']['vertical']['CEP95']:.4f}m/s")
            
            if 'attitude' in cep_stats:
                error_stats['attitude_cep'] = cep_stats['attitude']
                logging.info(f"姿态误差CEP: Pitch CEP50={cep_stats['attitude']['pitch']['CEP50']:.4f}°, Pitch CEP68={cep_stats['attitude']['pitch']['CEP68']:.4f}°, Pitch CEP95={cep_stats['attitude']['pitch']['CEP95']:.4f}°")
                logging.info(f"姿态误差CEP: Roll CEP50={cep_stats['attitude']['roll']['CEP50']:.4f}°, Roll CEP68={cep_stats['attitude']['roll']['CEP68']:.4f}°, Roll CEP95={cep_stats['attitude']['roll']['CEP95']:.4f}°")
                logging.info(f"姿态误差CEP: Yaw CEP50={cep_stats['attitude']['yaw']['CEP50']:.4f}°, Yaw CEP68={cep_stats['attitude']['yaw']['CEP68']:.4f}°, Yaw CEP95={cep_stats['attitude']['yaw']['CEP95']:.4f}°")
    
    # 构建保存路径：ref_file_path的相对路径
    ref_dir = os.path.dirname(ref_file_path)        
    save_dir = os.path.join(ref_dir, '..', '..', '..', 'batch_process', 'figure', data_name)
    save_dir = os.path.normpath(save_dir)  # 规范化路径
    
    # 创建保存目录
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        logging.info(f"✓ 创建保存目录: {save_dir}")
    
    # 保存所有图像到文件
    for plot_type, fig in figures:
        # 根据处理模式和图像类型生成文件名
        if 'LC' in process_mode or 'TC' in process_mode:
            # GNSS/INS模式：5张图片
            if plot_type == 'bg':
                filename = f"{data_name_without_ext}_bg.png"
            elif plot_type == 'ba':
                filename = f"{data_name_without_ext}_ba.png"
            elif plot_type == 'pos':
                filename = f"{data_name_without_ext}_pos.png"
            elif plot_type == 'vel':
                filename = f"{data_name_without_ext}_vel.png"
            elif plot_type == 'att':
                filename = f"{data_name_without_ext}_att.png"
        else:
            # GNSS模式：1张位置误差图
            filename = f"{data_name_without_ext}_pos.png"
        
        filepath = os.path.join(save_dir, filename)
        
        # 保存图像
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        logging.info(f"✓ 保存图像: {filepath}")
        
        # 关闭图形以释放内存
        plt.close(fig)
    
    logging.info(f"✓ 所有图像已保存到 {save_dir}")
    
    return error_stats