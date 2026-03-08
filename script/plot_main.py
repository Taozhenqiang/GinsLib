# main.py
import os
from read_sol import read_solution
from read_ref import read_ref
from plot_trajectory import plot_trajectory
from plot_position   import plot_position
from plot_velocity   import plot_velocity
from plot_attitude   import plot_attitude
from plot_imubias    import plot_imubias
from plot_err    import plot_err
import matplotlib.pyplot as plt

def main():
    # ================= 配置区域 =================
    # 在这里直接设置您要绘制的图表类型
    # 可用选项: 'trj'(轨迹), 'pos'(位置), 'vel'(速度), 'att'(姿态), 'bias'(零偏), 'err'(误差)
    PLOT_OPTIONS = ['trj','err']  # 修改这里来选择要绘制的图表
    
    # 误差类型配置（仅当选择err时使用）
    # 可用选项: 'p'(位置误差), 'v'(速度误差), 'a'(姿态误差)
    ERROR_TYPE = 'p'
    
    # 文件路径配置（可选，如果使用默认路径则保持为空）
    SOLFILE_PATH = ''  # 自定义结果文件路径，留空使用默认路径
    REFFILE_PATH = ''  # 自定义参考文件路径，留空使用默认路径

    # 图片保存配置
    SAVE_IMAGES = True  # 设置为True保存图像，False不保存
    # ================= 配置结束 =================
    
    # 检查配置是否有效
    valid_options = ['trj', 'pos_', 'vel_', 'att_', 'bias', 'err']
    for option in PLOT_OPTIONS:
        if option not in valid_options:
            print(f"错误: 无效的绘制选项 '{option}'")
            print(f"可用选项: {valid_options}")
            return
    
    if ERROR_TYPE not in ['p', 'v', 'a','pv','pa','va','pva']:
        print(f"错误: 无效的误差类型 '{ERROR_TYPE}'")
        print("可用选项: 'p'(位置误差), 'v'(速度误差), 'a'(姿态误差)")
        return

    # 提供文件路径和结果文件名
    if SOLFILE_PATH:
        solfile = SOLFILE_PATH
    else:
        pathname = 'F:/Navigation_Learn/GNSS_INS/GinsLib/data/GNSS_INS_Vehicle/STIM300_Vehicle_complex_20221230/result'
        filename = 'GNSS_Vehicle_complex_20221230_PPK_F.pos'
        solfile = os.path.join(pathname, filename)

    # 提供文件路径和参考文件名
    if REFFILE_PATH:
        refile = REFFILE_PATH
    else:
        pathname = 'F:/Navigation_Learn/GNSS_INS/GinsLib/data/GNSS_INS_Vehicle/STIM300_Vehicle_complex_20221230'
        filename = 'truth.truth'
        refile = os.path.join(pathname, filename)

    # 结果文件读取配置
    skip_lines = 28
    row = 10000
    # GINLIB: gps week, sow, pos[x/y/z], ratio, vel[x/y/z], att[pitch/roll/heading], bg[x/y/z], ba[x/y/z]
    col = 18
    sample = 1/100

    # 读取结果文件
    data = read_solution(solfile, skip_lines, row, col, sample)

    # 检查是否成功读取
    if data is not None:
        print("结果文件数据读取成功")
    else:
        print("结果文件数据读取失败")
        return

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
        elif plot_type == 'bias':
            print("绘制IMU零偏图...")
            fig_bg, fig_ba = plot_imubias(data)
            figures.append(('bg', fig_bg))
            figures.append(('ba', fig_ba))
        elif plot_type == 'err':
            print(f"绘制误差图 (类型: {ERROR_TYPE})...")
            fig, rms_stats = plot_err(data, ref_data, ERROR_TYPE)
            figures.extend(fig)  # 添加所有误差图
    
    # 保存所有图像到文件
    if SAVE_IMAGES:
        # 构建保存路径：参考文件同级的figure文件夹
        ref_dir = os.path.dirname(refile)
        save_dir = os.path.join(ref_dir, 'figure')
        
        # 创建保存目录
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            print(f"✓ 创建保存目录: {save_dir}")
        
        # 提取结果文件名（不含路径和扩展名）
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
            elif plot_type == 'att':
                filename = f"{data_name_without_ext}_att.png"
            elif plot_type == 'bg':
                filename = f"{data_name_without_ext}_bg.png"
            elif plot_type == 'ba':
                filename = f"{data_name_without_ext}_ba.png"
            elif plot_type == 'pos':
                filename = f"{data_name_without_ext}_pos_err.png"
            elif plot_type == 'vel':
                filename = f"{data_name_without_ext}_vel_err.png"
            elif plot_type == 'att':
                filename = f"{data_name_without_ext}_att_err.png"
            
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