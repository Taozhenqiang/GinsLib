import subprocess
import platform
import os
import time
import logging
import sys
from datetime import datetime
import pandas as pd

# 导入误差分析模块
sys.path.append(os.path.dirname(__file__))
from plot_batch import batch_plot_analysis

def setup_logging(data_root_path):
    """设置日志记录，同时输出到控制台和文件"""
    save_path = os.path.join(data_root_path, '..', '..','batch_process','log')
    save_path = os.path.normpath(save_path)
    
    # 创建保存目录
    if not os.path.exists(save_path):
        os.makedirs(save_path)
        print(f"✓ 创建保存目录: {save_path}")

    # 生成带时间戳的日志文件名
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(save_path , f"batch_process_{timestamp}.log")
    
    # 配置logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),  # 文件处理器
            logging.StreamHandler(sys.stdout)  # 控制台处理器
        ]
    )
    
    return log_file

def get_soltype(solt_value):
    """根据solt字段确定滤波模式"""
    if solt_value == 0:
        return "F"
    elif solt_value == 1:
        return "B"
    elif solt_value == 2:
        return "FB"
    else:
        return "F"  # 默认值

def get_result_filename(current_work_dir, mode, solt_value):
    """
    根据规则生成结果文件名
    
    参数:
    - current_work_dir: 当前工作目录路径
    - mode: 处理模式字符串
    - solt_value: solt字段值
    
    返回:
    - 结果文件名
    """
    # 提取文件夹名称
    folder_name = os.path.basename(current_work_dir)
    
    # 确定滤波模式
    soltype = get_soltype(solt_value)
    
    # 判断是否为GNSS/INS模式
    if 'LC' in mode or 'TC' in mode:
        # GNSS/INS模式：current_work_dir + "_" + mode中的第一部分 + "_" + soltype + 组合模式(LC或TC) + ".pos"
        mode_parts = mode.split()
        mode_first_part = mode_parts[0]  # 例如SPP、PPP等
        combined_mode = "LC" if 'LC' in mode else "TC"
        filename = f"{folder_name}_{mode_first_part}_{soltype}_{combined_mode}.pos"
    else:
        # GNSS模式："GNSS" + current_work_dir首个_符号后的部分 + "_" + mode + soltype + ".pos"
        # 提取首个_符号后的部分
        if '_' in folder_name:
            after_underscore = folder_name.split('_', 1)[1]
        else:
            after_underscore = folder_name
        filename = f"GNSS_{after_underscore}_{mode}_{soltype}.pos"
    
    return filename

def save_error_stats_to_excel(all_error_stats, save_path, mode):
    """
    将误差统计结果保存到Excel文件
    
    参数:
    - all_error_stats: 所有误差统计结果的列表
    - save_path: 保存路径
    """
    # 创建数据列表
    data_rows = []
    
    for stats in all_error_stats:
        # 创建一行数据
        row_data = {
            '数据': stats['data_name'],
            '模式': stats['process_mode']
        }
        
        # 位置误差
        if stats['position_rms']:
            pos_rms = stats['position_rms']
            row_data['E (m)'] = pos_rms['E']
            row_data['N (m)'] = pos_rms['N']
            row_data['U (m)'] = pos_rms['U']
            row_data['3D (m)'] = pos_rms['3D']
        else:
            row_data['E (m)'] = None
            row_data['N (m)'] = None
            row_data['U (m)'] = None
            row_data['3D (m)'] = None
        
        # 速度误差（仅GNSS/INS模式）
        if 'LC' in stats['process_mode'] or 'TC' in stats['process_mode']:
            if stats['velocity_rms']:
                vel_rms = stats['velocity_rms']
                row_data['E (m/s)'] = vel_rms['E']
                row_data['N (m/s)'] = vel_rms['N']
                row_data['U (m/s)'] = vel_rms['U']
                row_data['3D (m/s)'] = vel_rms['3D']
            else:
                row_data['E (m/s)'] = None
                row_data['N (m/s)'] = None
                row_data['U (m/s)'] = None
                row_data['3D (m/s)'] = None
            
            # 姿态误差（仅GNSS/INS模式）
            if stats['attitude_rms']:
                att_rms = stats['attitude_rms']
                row_data['pitch (deg)'] = att_rms['pitch']
                row_data['roll (deg)'] = att_rms['roll']
                row_data['yaw (deg)'] = att_rms['yaw']
            else:
                row_data['pitch (deg)'] = None
                row_data['roll (deg)'] = None
                row_data['yaw (deg)'] = None
        else:
            # GNSS模式，速度误差和姿态误差为空
            row_data['E (m/s)'] = None
            row_data['N (m/s)'] = None
            row_data['U (m/s)'] = None
            row_data['3D (m/s)'] = None
            row_data['pitch (deg)'] = None
            row_data['roll (deg)'] = None
            row_data['yaw (deg)'] = None
        
        data_rows.append(row_data)
    
    # 创建DataFrame
    df = pd.DataFrame(data_rows)
    
    # 定义列的顺序
    columns_order = [
        '数据', '模式', 
        'E (m)', 'N (m)', 'U (m)', '3D (m)',
        'E (m/s)', 'N (m/s)', 'U (m/s)', '3D (m/s)',
        'pitch (deg)', 'roll (deg)', 'yaw (deg)'
    ]
    
    # 重新排列列的顺序
    df = df.reindex(columns=columns_order)
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    # 生成Excel文件名（带时间戳）
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    excel_file = os.path.join(save_path, f"error_stats_{mode}_{timestamp}.xlsx")
    
    # 保存到Excel文件
    df.to_excel(excel_file, index=False, engine='openpyxl')
    
    return excel_file

def batch_process_gins(exe_path_input):

    # 1. 获取当前工作目录 (对应 VSCode launch.json 中的 "cwd")
    # 假设当前 cwd 就是 data_gins 文件夹
    data_root_path = os.getcwd()

    # 设置日志记录
    log_file = setup_logging(data_root_path)
    logging.info(f"日志文件保存位置: {log_file}")
    
    # 记录批处理开始时间
    start_time = time.time()
    logging.info("批处理任务开始执行")
    
    # 将可执行文件路径转换为绝对路径
    # 注意：这里的 exe_path_input 是相对于 cwd 的路径
    exe_path = os.path.abspath(exe_path_input)
    
    logging.info(f"当前工作目录 (Data Root): {data_root_path}")
    logging.info(f"解算程序路径 (EXE Path):  {exe_path}")
    logging.info("-" * 60)

    # 检查可执行文件是否存在
    if not os.path.exists(exe_path):
        logging.error(f"无法找到可执行文件，请检查路径配置: {exe_path}")
        return

    # 2. 遍历当前目录下的所有子文件夹
    # 过滤掉非文件夹项
    sub_folders = [f for f in os.listdir(data_root_path) if os.path.isdir(os.path.join(data_root_path, f))]
    
    if not sub_folders:
        logging.warning("当前目录下没有发现子文件夹。")
        return

    logging.info(f"共 {len(sub_folders)} 个数据任务，开始批处理...")

    # 存储所有误差统计结果
    all_error_stats = []

    # for folder_name in sub_folders:  
    for folder_name in sub_folders[31:32]:  # 处理单个子文件夹用于测试
        # 进入子文件夹路径
        current_work_dir = os.path.join(data_root_path, folder_name)
        
        # 3. 【修改点】配置文件名默认与数据文件夹名一致
        # 例如文件夹为 "data1"，则寻找 "data1.conf"
        conf_file = f"{folder_name}.conf"
        target_conf_path = os.path.join(current_work_dir, conf_file)
        
        logging.info(f"第 {sub_folders.index(folder_name) + 1} 个任务: {folder_name}")
        logging.info(f"- 路径: {current_work_dir}")
        logging.info(f"- 配置: {conf_file}")

        # 检查配置文件是否存在
        if not os.path.exists(target_conf_path):
            logging.error(f"配置文件不存在: {target_conf_path}")
            logging.info("跳过该任务。")
            continue
        
        # 4. 定义处理模式 (保持原有逻辑)
        # mode = ["SPP","SPP LC","SPP TC","PPP","PPP LC","PPP TC","PPD","PPD LC","PPD TC","PPK","PPK LC","PPK TC"]
        mode = ["SPP","SPP LC","SPP TC","PPD","PPD LC","PPD TC"]
        
        # 命令列表
        commands = [
            [exe_path,"-k",conf_file,"-p","0","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0","-solt","0"], #spp
            [exe_path,"-k",conf_file,"-p","0","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0","-solt","0"], #spp/ins LC
            [exe_path,"-k",conf_file,"-p","0","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0","-solt","0"], #spp/ins TC
            # [exe_path,"-k",conf_file,"-p","8","-gins","0","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0","-solt","0"], #ppp
            # [exe_path,"-k",conf_file,"-p","8","-gins","1","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0","-solt","0"], #ppp/ins LC
            # [exe_path,"-k",conf_file,"-p","8","-gins","2","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0","-solt","0"], #ppp/ins TC
            [exe_path,"-k",conf_file,"-p","1","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0","-solt","0"], #ppd
            [exe_path,"-k",conf_file,"-p","1","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0","-solt","0"], #ppd/ins LC
            [exe_path,"-k",conf_file,"-p","1","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0","-solt","0"], #ppd/ins TC
            # [exe_path,"-k",conf_file,"-p","2","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2","-solt","0"], #ppk
            # [exe_path,"-k",conf_file,"-p","2","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2","-solt","0"], #ppp/ins LC
            # [exe_path,"-k",conf_file,"-p","2","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2","-solt","0"], #ppp/ins TC
        ]

        # 5. 循环执行命令
        all_success = True
        for command, m in zip(commands, mode):
            try:
                # 提取solt值用于文件名生成
                solt_index = command.index("-solt") + 1
                solt_value = int(command[solt_index])
                
                # cwd=current_work_dir 确保程序在子文件夹内部运行，从而能读取到 conf_file   
                process = subprocess.Popen(command, cwd=current_work_dir, 
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                          text=True, bufsize=1, universal_newlines=True)
                
                # 实时读取并打印输出，添加去重机制
                last_line = None
                for line in process.stderr:
                    stripped_line = line.strip()
                    # 跳过空行和与上一行相同的行
                    if stripped_line and stripped_line != last_line:
                        logging.info(f"[GINSLIB] {stripped_line}")
                        last_line = stripped_line
                
                # 等待进程结束
                process.wait()
                result = process
                
                if result.returncode == 0:
                    logging.info(f"✓ {folder_name} {m} 执行成功")
                    
                    # 6. 执行误差分析
                    try:
                        # 确定结果文件路径
                        result_filename = get_result_filename(current_work_dir, m, solt_value)
                        result_file_path = os.path.join(current_work_dir, "result", result_filename)
                        
                        # 确定参考文件路径
                        ref_file_path = os.path.join(current_work_dir, "truth.truth")
                        
                        # 检查文件是否存在
                        if not os.path.exists(result_file_path):
                            logging.warning(f"结果文件不存在，跳过误差分析: {result_file_path}")
                            continue
                            
                        if not os.path.exists(ref_file_path):
                            logging.warning(f"参考文件不存在，跳过误差分析: {ref_file_path}")
                            continue
                        
                        logging.info(f"开始误差分析")
                        logging.info(f"- 结果文件: {result_file_path}")
                        logging.info(f"- 参考文件: {ref_file_path}")
                        
                        # 调用误差分析函数
                        error_stats = batch_plot_analysis(result_file_path, ref_file_path, m)
                        all_error_stats.append(error_stats)
                        
                        logging.info(f"✓ 误差分析完成: {result_filename}")
                        
                    except Exception as e:
                        logging.error(f"误差分析异常: {e}")
                        
                else:
                    all_success = False
                    logging.error(f"✗ {m} 执行失败 (返回码: {result.returncode})")
                    # 如有报错，输出部分错误信息
                    if result.stderr:
                         logging.error(f"错误信息: {result.stderr.strip()[:200]}")
            except Exception as e:
                all_success = False
                logging.error(f"执行异常: {e}")

    # 保存误差统计结果到Excel文件
    if all_error_stats:
        # 确定保存路径：./current_work_dir../batch_process
        save_path = os.path.join(data_root_path, '..', '..', 'batch_process', 'excel')
        save_path = os.path.normpath(save_path)
        
        try:
            excel_file = save_error_stats_to_excel(all_error_stats, save_path, m)
            logging.info(f"✓ 误差统计结果保存到Excel文件: {excel_file}")
        except Exception as e:
            logging.error(f"保存Excel文件失败: {e}")
    else:
        logging.warning("没有误差统计结果需要保存")

    # 计算总耗时
    end_time = time.time()
    total_time = end_time - start_time
    
    logging.info(f"{len(sub_folders)} 个批处理任务结束。")
    logging.info(f"总共耗费时间: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")
    logging.info(f"日志文件保存路径: {log_file}")

if __name__ == "__main__":
    # ================= 配置区域 =================
    
    # 设置 GINSLIB 可执行文件的路径
    # 请注意：这里的路径是相对于 data_gins (即 cwd) 的
    # 相对路径，或者是绝对路径
    # 假设你的目录结构是:
    # workspace/
    #   ├── build/Bin/GINSLIB.exe
    #   └── data_gins/  <-- VSCode cwd 在这里
    #         ├── data1/
    #         └── ...
    # 那么相对路径应该是 "../build/Bin/GINSLIB.exe"
    
    EXE_PATH = r"../../build/Bin/GINSLIB.exe" 
    
    # 操作系统适配
    if platform.system() != "Windows":
        EXE_PATH = "../build/Bin/GINSLIB"

    # ================= 开始运行 =================
    batch_process_gins(EXE_PATH)