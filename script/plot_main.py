# main.py
import os
from read_file import read_pos
from read_file import read_ref
from plot_solution import plot_trajectory
from plot_solution import plot_position
from plot_solution import plot_velocity
from plot_solution import plot_attitude
from plot_solution import plot_nsatdop
from plot_solution import plot_ratio
from plot_solution import plot_solflag
from plot_solution import plot_pvastd
from plot_solution import plot_imubias
from plot_err    import plot_err, plot_err_multi
import matplotlib.pyplot as plt

def main():
    # ================= 配置区域 =================
    # 在这里直接设置您要绘制的图表类型
    # 可用选项: 'trj'(轨迹), 'pos_(位置), 'vel_'(速度), 'att_(姿态), 'nsat(卫星数)', 'ratio', 'solflag', 'pvastd', 'bias'(零偏), 'err'(误差)
    PLOT_OPTIONS = ['solflag','err']  # 修改这里来选择要绘制的图表
    
    # 误差类型配置（仅当选择err时使用）
    # 可用选项: 'p'(位置误差), 'v'(速度误差), 'a'(姿态误差)
    ERROR_TYPE = 'p'
    
    # 文件路径配置
    PATH_NAME = './GNSS/LG69T_Vehicle_complex_20250414'
    
    # 多文件误差分析选项
    MULTI_FILE_ANALYSIS = False  # True: 多文件对比分析, False: 单文件分析
    
    # 单文件分析配置
    SOLFILE_PATH = './result/GNSS_Vehicle_complex_20250414_SPP_F.pos'  # pos文件路径
    
    # 多文件分析配置（当MULTI_FILE_ANALYSIS为True时使用）
    SOLFILE_PATHS = [
        './result/GNSS_Vehicle_complex_20250414_PPK_F_INST.pos',
        './result/ASM330_Vehicle_complex_20250414_GNSS_F_LC.pos', 
    ]
    SOLFILE_LABELS = [
        'PPK',
        'PPK/INS-LC',
    ]
    TIME_SYNC_MODE = 1  # 时间同步模式 (1: 各自独立时间序列; 2: 以第一个文件为基准同步)

    REFFILE_PATH = 'truth_pva.truth'  # ref参考文件路径

    # pos文件配置
    SOL_TYPE = 'GINSLIB'  # pos文件类型，GINSLIB/OTHER
    SOL_IDX  = [0,1,2,3,4,5,6,7,9,8,10]         # 第三方结果文件列索引
    SOLPOS_TYPE = 'LLH'  # pos文件坐标类型，ECEF/LLH
    SOLVEL_TYPE = 'NED'  # vel文件坐标类型，ECEF/ENU/NED

    # 图片保存配置
    SAVE_IMAGES = False  # 设置为True保存图像，False不保存
    DEFINE_PATH = True
    SAVE_PATH =  './figure/PPK/GECJ+CONT_FLOAT'
    # ================= 配置结束 =================
    
    # 检查配置是否有效
    valid_options = ['trj', 'pos_', 'vel_', 'att_', 'nsat', 'ratio', 'solflag', 'pvastd', 'bias', 'err']
    for option in PLOT_OPTIONS:
        if option not in valid_options:
            print(f"错误: 无效的绘制选项 '{option}'")
            print(f"可用选项: {valid_options}")
            return
    
    if ERROR_TYPE not in ['p', 'v', 'a','pv','pa','va','pva']:
        print(f"错误: 无效的误差类型 '{ERROR_TYPE}'")
        print("可用选项: 'p'(位置误差), 'v'(速度误差), 'a'(姿态误差)")
        return
    
    # 检查多文件分析配置
    if MULTI_FILE_ANALYSIS:
        if len(SOLFILE_PATHS) != len(SOLFILE_LABELS):
            print("错误: 多文件分析时，SOLFILE_PATHS和SOLFILE_LABELS数量必须相同")
            return
        if len(SOLFILE_PATHS) == 0:
            print("错误: 多文件分析时，必须提供至少一个结果文件路径")
            return

    pathname = PATH_NAME
    
    # 多文件分析模式
    if MULTI_FILE_ANALYSIS:
        # 读取多个结果文件
        solutions = []
        for sol_path in SOLFILE_PATHS:
            solfile = os.path.join(pathname, sol_path)
            # 结果文件读取配置
            skip_lines = 0
            row = 10000
            # GINLIB: gps week, sow, pos[x/y/z], ratio, vel[x/y/z], att[pitch/roll/heading], bg[x/y/z], ba[x/y/z]
            col = 36
            sample = 1/100
            data = read_pos(solfile, skip_lines, row, col, sample, SOL_TYPE, SOL_IDX, SOLPOS_TYPE, SOLVEL_TYPE)
            if data is not None:
                solutions.append(data)
                print(f"结果文件 {sol_path} 读取成功")
            else:
                print(f"结果文件 {sol_path} 读取失败")
                solutions.append(None)
        
        # 过滤掉读取失败的文件
        valid_solutions = []
        valid_labels = []
        for i, (sol, label) in enumerate(zip(solutions, SOLFILE_LABELS)):
            if sol is not None:
                valid_solutions.append(sol)
                valid_labels.append(label)
        
        if len(valid_solutions) == 0:
            print("错误: 所有结果文件读取失败")
            return
            
        print(f"成功读取 {len(valid_solutions)} 个结果文件")
    else:
        # 单文件分析模式
        if SOLFILE_PATH:
            solfile = os.path.join(pathname, SOLFILE_PATH)
        else:
            filename = './result/ASM330_Vehicle_complex_20250414_SPP_F_LC.pos'
            solfile = os.path.join(pathname, filename)

        # 结果文件读取配置
        skip_lines = 0
        row = 10000
        # GINLIB: gps week, sow, pos[x/y/z], ratio, vel[x/y/z], att[pitch/roll/heading], bg[x/y/z], ba[x/y/z]
        col = 36
        sample = 1/100

        # 读取结果文件
        data = read_pos(solfile, skip_lines, row, col, sample, SOL_TYPE, SOL_IDX, SOLPOS_TYPE, SOLVEL_TYPE)

        # 检查是否成功读取
        if data is not None:
            print("结果文件数据读取成功")
        else:
            print("结果文件数据读取失败")
            return

    # 提供文件路径和参考文件名
    if REFFILE_PATH:
        refile = os.path.join(pathname, REFFILE_PATH)
    else:
        filename = 'truth_pva.truth'
        refile = os.path.join(pathname, filename)

    # 参考文件读取配置
    skip_lines = 20
    row = 10000
    # gps week, sow, pos[x/y/z], vel[x/y/z], att[pitch/roll/heading]
    col = 11
    ins_flag = 1
    ins_interval = 1/100
    GNSS_interval = 1
    # truth-tzq
    idx = [0,1,2,3,4,5,6,7,8,9,10]
    # IE-tzq
    # idx = [0,1,2,3,4,8,9,10,14,15,16]
    # IE-zf
    # idx = [0,1,9,10,11,15,16,17,22,23,21]

    if 'err' in PLOT_OPTIONS:
        # 读取参考文件
        ref_data = read_ref(refile, skip_lines, row, col, ins_flag, ins_interval, GNSS_interval, idx)

        # 检查是否成功读取
        if ref_data is not None:
            print("参考文件数据读取成功")
        else:
            print("参考文件数据读取失败")
            return
    
    
    # 根据配置选择绘制图表
    print(f"开始绘制选定的图表: {PLOT_OPTIONS}")
    
    # 在绘制循环中修改
    figures = []  # 存储所有图形对象和对应的类型

    for plot_type in PLOT_OPTIONS:
        if plot_type == 'trj':
            print("绘制轨迹图...")
            fig = plot_trajectory(data)
            figures.append(('trj', fig))
        elif plot_type == 'pos_':
            print("绘制位置图...")
            fig = plot_position(data)
            figures.append(('pos_', fig))
        elif plot_type == 'vel_':
            print("绘制速度图...")
            fig = plot_velocity(data)
            figures.append(('vel_', fig))
        elif plot_type == 'att_':
            print("绘制姿态图...")
            fig = plot_attitude(data)
            figures.append(('att_', fig))
        elif plot_type== 'nsat':
            fig = plot_nsatdop(data)
            figures.append(('nsat', fig))
        elif plot_type == 'ratio':
            fig = plot_ratio(data)
            figures.append(('ratio', fig))
        elif plot_type == 'solflag':
            fig = plot_solflag(data)
            figures.append(('solflag', fig))
        elif plot_type == 'pvastd':
            fig = plot_pvastd(data)
            figures.append(('pvastd', fig))
        elif plot_type == 'bias':
            print("绘制IMU零偏图...")
            fig_bg, fig_ba = plot_imubias(data)
            figures.append(('bg', fig_bg))
            figures.append(('ba', fig_ba))
        elif plot_type == 'err':
            print(f"绘制误差图 (类型: {ERROR_TYPE})...")
            
            if MULTI_FILE_ANALYSIS:
                # 多文件误差分析
                print("使用多文件误差分析模式...")
                fig, rms_stats, cep_stats = plot_err_multi(valid_solutions, valid_labels, ref_data, ERROR_TYPE, TIME_SYNC_MODE)
                figures.extend(fig)  # 添加所有误差图
            else:
                # 单文件误差分析
                fig, rms_stats, cep_stats = plot_err(data, ref_data, ERROR_TYPE, True)
                figures.extend(fig)  # 添加所有误差图
    
    # 保存所有图像到文件
    if SAVE_IMAGES:
        if DEFINE_PATH:
            # 构建保存路径：参考文件同级的figure文件夹
            save_dir = os.path.join(pathname, SAVE_PATH)
        else:
            save_dir = os.path.join(pathname, 'figure')
        
        # 创建保存目录
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            print(f"✓ 创建保存目录: {save_dir}")
        
        # 提取结果文件名（不含路径和扩展名）
        if MULTI_FILE_ANALYSIS:
            # 多文件分析使用第一个文件的名称
            first_sol_path = SOLFILE_PATHS[0]
            sol_filename = os.path.basename(first_sol_path)
        else:
            sol_filename = os.path.basename(solfile)
            
        if sol_filename.endswith('.pos'):
            data_name_without_ext = sol_filename[:-4]  # 去掉.pos扩展名
        else:
            data_name_without_ext = sol_filename
        
        # 保存所有图像到文件
        for plot_type, fig in figures:
            # 根据图像类型生成文件名
            if plot_type == 'trj':
                filename = f"{data_name_without_ext}_trj.png"
            elif plot_type == 'pos_':
                filename = f"{data_name_without_ext}_pos.png"
            elif plot_type == 'vel_':
                filename = f"{data_name_without_ext}_vel.png"
            elif plot_type == 'att_':
                filename = f"{data_name_without_ext}_att.png"
            elif plot_type == 'bg':
                filename = f"{data_name_without_ext}_bg.png"
            elif plot_type == 'ba':
                filename = f"{data_name_without_ext}_ba.png"
            elif plot_type == 'pos':
                filename = f"{data_name_without_ext}_pos_err.png"
            elif plot_type == '3Dpos':
                filename = f"{data_name_without_ext}_3Dpos_err.png"
            elif plot_type == 'pos_cdf':
                filename = f"{data_name_without_ext}_pos_cdf.png"
            elif plot_type == 'vel':
                filename = f"{data_name_without_ext}_vel_err.png"
            elif plot_type == 'att':
                filename = f"{data_name_without_ext}_att_err.png"
            elif plot_type == 'pos_multi':
                filename = f"{data_name_without_ext}_multi_pos_err.png"
            elif plot_type == 'vel_multi':
                filename = f"{data_name_without_ext}_multi_vel_err.png"
            elif plot_type == 'att_multi':
                filename = f"{data_name_without_ext}_multi_att_err.png"
            
            filepath = os.path.join(save_dir, filename)
            
            # 保存图像
            fig.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"✓ 保存图像: {filepath}")
        
        print(f"✓ 所有图像已保存到 {save_dir}")
    else:
        print("图片保存已关闭，仅显示图像")
    
    # 最后统一显示所有图形
    print("所有图表绘制完成，开始显示...")
    for plot_type, fig in figures:
        plt.show(block=False)   # 非阻塞显示

    # 保持程序运行，直到所有窗口关闭
    plt.show(block=True)

if __name__ == "__main__":
    main()