import subprocess
import platform
import os
import time
import logging
import sys
from datetime import datetime

def setup_logging():
    """设置日志记录，同时输出到控制台和文件"""
    data_root_path = os.getcwd()
    
    # 生成带时间戳的日志文件名
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(data_root_path, f"batch_process_{timestamp}.log")
    
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

def batch_process_gins(exe_path_input):
    # 设置日志记录
    log_file = setup_logging()
    logging.info(f"日志文件保存位置: {log_file}")
    
    # 记录批处理开始时间
    start_time = time.time()
    logging.info("批处理任务开始执行")
    
    # 1. 获取当前工作目录 (对应 VSCode launch.json 中的 "cwd")
    # 假设当前 cwd 就是 data_gins 文件夹
    data_root_path = os.getcwd()
    
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

    for folder_name in sub_folders:
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
        mode = ["spp","spp/ins LC","spp/ins TC","ppp","ppp/ins LC","ppp/ins TC","ppd","ppd/ins LC","ppd/ins TC","ppk","ppk/ins LC","ppk/ins TC"]
        
        # 命令列表
        commands = [
            [exe_path,"-k",conf_file,"-p","0","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #spp
            # [exe_path,"-k",conf_file,"-p","0","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #spp/ins LC
            # [exe_path,"-k",conf_file,"-p","0","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #spp/ins TC
            # [exe_path,"-k",conf_file,"-p","8","-gins","0","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0"], #ppp
            # [exe_path,"-k",conf_file,"-p","8","-gins","1","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0"], #ppp/ins LC
            # [exe_path,"-k",conf_file,"-p","8","-gins","2","-ion","4","-tro","3","-eph","1","-flt","0","-amb","0"], #ppp/ins TC
            # [exe_path,"-k",conf_file,"-p","1","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #ppd
            # [exe_path,"-k",conf_file,"-p","1","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #ppp/ins LC
            # [exe_path,"-k",conf_file,"-p","1","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","0"], #ppp/ins TC
            # [exe_path,"-k",conf_file,"-p","2","-gins","0","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2"], #ppk
            # [exe_path,"-k",conf_file,"-p","2","-gins","1","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2"], #ppp/ins LC
            # [exe_path,"-k",conf_file,"-p","2","-gins","2","-ion","1","-tro","1","-eph","0","-flt","0","-amb","2"], #ppp/ins TC
        ]

        # 5. 循环执行命令
        all_success = True
        for command, m in zip(commands, mode):
            try:
                # cwd=current_work_dir 确保程序在子文件夹内部运行，从而能读取到 conf_file   
                process = subprocess.Popen(command, cwd=current_work_dir, 
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                          text=True, bufsize=1, universal_newlines=True)
                
                # 实时读取并打印输出
                for line in process.stderr:
                    # 同时记录到日志文件
                    logging.info(f"[GINSLIB] {line.strip()}")
                
                # 等待进程结束
                process.wait()
                result = process
                
                if result.returncode == 0:
                    logging.info(f"✓ {folder_name} {m} 执行成功")
                else:
                    all_success = False
                    logging.error(f"✗ {m} 执行失败 (返回码: {result.returncode})")
                    # 如有报错，输出部分错误信息
                    if result.stderr:
                         logging.error(f"错误信息: {result.stderr.strip()[:200]}")
            except Exception as e:
                all_success = False
                logging.error(f"执行异常: {e}")

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
    
    EXE_PATH = r"..\..\build\Bin\GINSLIB.exe" 
    
    # 操作系统适配
    if platform.system() != "Windows":
        EXE_PATH = "../build/Bin/GINSLIB"

    # ================= 开始运行 =================
    batch_process_gins(EXE_PATH)