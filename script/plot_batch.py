import os
import numpy as np
import logging
from read_sol import read_solution
from read_ref import read_ref
from plot_trajectory import plot_trajectory
from plot_position   import plot_position
from plot_velocity   import plot_velocity
from plot_attitude   import plot_attitude
from plot_imubias    import plot_imubias
from plot_err    import plot_err
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
    
    # 提取data_name：sol_file_path中倒数第三个分隔符与倒数第二个分隔符之间的部分
    path_parts = sol_file_path.split(os.sep)
    if len(path_parts) >= 3:
        # 获取倒数第三个分隔符与倒数第二个分隔符之间的部分
        data_name = path_parts[-3] if len(path_parts) >= 3 else path_parts[-1]
    else:
        # 如果路径分隔符不足3个，则使用文件名
        data_name = os.path.basename(sol_file_path)

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
        'attitude_rms': None
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
            err_figures, rms_stats = plot_err(data, ref_data, error_type)
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