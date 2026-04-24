import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from pathlib import Path

def read_imu_data(file_path, time_col, gyro_cols, accel_cols, time_format, skip_rows=0, end_row=0):
    """
    读取IMU数据文件
    args:
        file_path: 输入文件路径
        time_col: 时间列索引
        time_format: 时间格式：'s'(秒) 或 'ms'(毫秒)
        gyro_cols: 陀螺仪数据列索引列表 [x, y, z]
        accel_cols: 加速度计数据列索引列表 [x, y, z]
        skip_rows: 跳过的行数
        end_row: 结束行数（可选，留空读取到文件末尾）
    return: 时间戳数组，陀螺仪数据数组，加速度计数据数组
    """
    try:
        # 根据文件扩展名选择读取方法
        file_ext = Path(file_path).suffix.lower()

        if time_format.lower() == 'ms':
            sample_interval = 0.001
        elif time_format.lower() == 's':
            sample_interval = 1.0
        else:
            raise ValueError("时间格式必须为's'或'ms'")
        
        if file_ext in ['.csv']:
            df = pd.read_csv(file_path, skiprows=skip_rows, nrows=end_row-skip_rows+1, header=None, index=None)
        elif file_ext in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path, skiprows=skip_rows, nrows=end_row-skip_rows+1, header=None, index=None)
        else:
            # 文本文件处理
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                if end_row > 0:
                    lines = lines[:end_row]
                else:
                    lines = lines[skip_rows:]
            
            # 跳过指定行数
            data_lines = lines[skip_rows:end_row]
            
            # 解析数据
            data = []
            for line in data_lines:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('%'):
                    continue
                
                # 分割数据（支持空格、制表符、逗号分隔）
                elements = line.replace('\t', ' ').replace(',', ' ').split()
                try:
                    row_data = [float(elem) for elem in elements]
                    data.append(row_data)
                except ValueError:
                    continue
            
            if not data:
                raise ValueError("无法解析文件数据")
            
            df = pd.DataFrame(data)
        
        # 提取时间数据
        time_data = (df.iloc[:, time_col].values - df.iloc[0, time_col]) * sample_interval  # 以第一个时间戳为起点，转换为相对时间，并将时间尺度转换为秒
        
        # 提取陀螺仪数据
        gyro_data = df.iloc[:, gyro_cols].values
        
        # 提取加速度计数据
        accel_data = df.iloc[:, accel_cols].values
        
        return time_data, gyro_data, accel_data
        
    except Exception as e:
        print(f"读取文件错误: {e}")
        return None, None, None

def convert_gyro_units2rad(gyro_data, input_unit, data_format, sample_interval):
    """
    转换陀螺仪数据单位
    args:
        gyro_data: 陀螺仪原始数据
        input_unit: 输入单位 ('deg', 'rad', 'deg/s', 'rad/s')
        data_format: 数据格式 ('rate' 或 'increment')
        sample_interval: 采样间隔（秒）
    return: 转换后的陀螺仪数据（单位：rad）
    """
    converted_data = gyro_data.copy()
    
    # 单位转换
    if input_unit == 'deg':
        # 度转换为度/秒（如果已经是速率形式，则不需要转换）
        if data_format == 'increment':
            converted_data = np.deg2rad(converted_data)
    
    elif input_unit == 'rad':
        # 已经是目标单位，不需要转换
        pass
    
    elif input_unit == 'rad/s':
        # 弧度/秒转换为度/秒
        converted_data = converted_data * sample_interval  # rad/s * s = rad
    
    elif input_unit == 'deg/s':
        converted_data = np.deg2rad(converted_data) * sample_interval  # deg/s * s = deg → rad
    
    else:
        print(f"警告: 未知的陀螺仪单位 '{input_unit}'，假设为 deg/s")
    
    return converted_data

def convert_accel_units2mg(accel_data, input_unit):
    """
    转换加速度计数据单位
    args:
        accel_data: 加速度计原始数据
        input_unit: 输入单位 ('m/s^2' 或 'g')
    return: 转换后的加速度计数据（单位：mg）
    """
    converted_data = accel_data.copy()
    
    if input_unit == 'm/s^2':
        # m/s^2 转换为 g (1 g = 9.80665 m/s^2)
        converted_data = converted_data / 9.80665 * 1000  # 转换为 mg
    
    elif input_unit == 'g':
        # 已经是目标单位，不需要转换
        converted_data = converted_data * 1000  # 转换为 mg
    
    else:
        print(f"警告: 未知的加速度计单位 '{input_unit}'，假设为 g")
    
    return converted_data

def avar(y0, tau0, str=''):
    """
    Calculate Allan variance 
    输入：
        y0 - 输入数据 (陀螺: deg/h; 加速度计: mg)，一维数组/列表
        tau0 - 采样间隔 (s)
        str - 纵轴标签字符串
    输出：
        sigma - Allan方差值
        tau - 相关时间 (s)
        Err - Allan方差误差边界
    """
    
    # 初始化参数
    y0 = np.array(y0).flatten()  # 转为一维numpy数组
    N = len(y0)
    y = y0.copy()
    NL = N
    sigma = []
    tau = []
    Err = []
    
    # 核心Allan方差计算循环
    max_k = int(np.log2(N))
    for k in range(max_k):
        # 计算Allan标准差 
        diff = y[1:NL] - y[0:NL-1]
        sigma_k = np.sqrt(1 / (2 * (NL - 1)) * np.sum(diff ** 2))
        sigma.append(sigma_k)
        
        # 计算相关时间 tau
        tau_k = (2 ** k) * tau0
        tau.append(tau_k)
        
        # 计算误差边界
        err_k = 1 / np.sqrt(2 * (NL - 1))
        Err.append(err_k)
        
        # 数据长度减半 + 相邻点平均
        NL = NL // 2
        if NL < 3:
            break
        y = 0.5 * (y[0:2*NL:2] + y[1:2*NL:2])
    
    # 转为numpy数组
    sigma = np.array(sigma).reshape(-1, 1)
    tau = np.array(tau).reshape(-1, 1)
    Err = np.array(Err).reshape(-1, 1)
    
    return sigma, tau, Err


def avarimu(imu, sample_interval=0.01):
    """
    计算SIMU陀螺和加速度计的Allan方差
    输入：
        imu - IMU数据数组，格式：[gx,gy,gz,ax,ay,az, t] 
              列：1-3=陀螺角增量，4-6=加速度计速度增量，最后一列=时间戳
    输出：
        sigma - Allan方差 (6列：陀螺3轴+加速度计3轴)
        tau - 相关时间 (同列数)
    """
    dph = np.pi/180/3600    # deg/h → rad/s 转换系数
    mg  = 1e-3 * 9.80665     # mg → m/s² 转换系数

    # 初始化变量
    imu = np.array(imu, dtype=np.float64)
    sigma = []
    tau = []
    m = np.zeros(6)  # 存储陀螺、加速度计均值
    
    # 计算采样间隔 ts (对应MATLAB diff(imu(1:2,end)))
    ts = sample_interval  # 直接使用配置的采样间隔
    
    # 处理陀螺3轴 (MATLAB k=1:3 → Python 0:2)
    for k in range(3):
        g = imu[:, k] / ts / dph  # 单位转换：角增量 → deg/h
        m[k] = np.mean(g)             # 计算均值
        imu[:, k] = g - m[k]          # 去均值
        # 调用avar计算，关闭绘图
        sig_k, tau_k, _ = avar(g, ts, str='')
        sigma.append(sig_k.flatten())
        tau.append(tau_k.flatten())
    
    # 处理加速度计3轴 (MATLAB k=4:6 → Python 3:5)
    for k in range(3, 6):
        a = imu[:, k] / ts / mg   # 单位转换：速度增量 → mg
        m[k] = np.mean(a)             # 计算均值
        imu[:, k] = a - m[k]          # 去均值
        # 调用avar计算，关闭绘图
        sig_k, tau_k, _ = avar(a, ts, str='')
        sigma.append(sig_k.flatten())
        tau.append(tau_k.flatten())
    
    # 转置为列向量格式
    sigma = np.array(sigma).T
    tau = np.array(tau).T

    # 绘图
    plt.figure(figsize=(10, 8), facecolor='white')
    
    # 子图1：陀螺原始数据 (左上)
    plt.subplot(2, 2, 1)
    plt.plot(imu[:, -1], imu[:, 0:3])
    plt.grid(True)
    plt.xlabel(r'$t$ / s')
    plt.ylabel(r'$\omega$ / (°)/h')

    # 子图2：加速度计原始数据 (右上)
    plt.subplot(2, 2, 2)
    plt.plot(imu[:, -1], imu[:, 3:6])
    plt.grid(True)
    plt.xlabel(r'$t$ / s')
    plt.ylabel(r'$f^b$ / mg')

    # 子图3：陀螺Allan方差 (左下，双对数坐标)
    plt.subplot(2, 2, 3)
    plt.loglog(tau[:, 0], sigma[:, 0:3])
    plt.grid(True, which="both", ls="-")
    plt.xlabel(r'$\tau$ / s')
    plt.ylabel(r'$\sigma_A(\tau)$ / (°)/h')

    # 子图4：加速度计Allan方差 (右下，双对数坐标)
    plt.subplot(2, 2, 4)
    plt.loglog(tau[:, 0], sigma[:, 3:6])
    plt.grid(True, which="both", ls="-")
    plt.xlabel(r'$\tau$ / s')
    plt.ylabel(r'$\sigma_A(\tau)$ / mg')

    plt.tight_layout()
    
    return sigma, tau

def main():
    """主函数"""
    # ================= 配置区域 =================
    # 在这里直接设置IMU数据文件的参数配置
    
    # 文件路径配置（必需）
    INPUT_FILE = './GNSS/LUA300C/LUA300C(00).txt'  # 输入文件路径
    
    # 数据列索引配置（必需）
    TIME_COL = 0                        # 时间列索引（从0开始）
    GYRO_COLS = [2, 3, 4]               # 陀螺仪数据列索引 [x, y, z]
    ACCEL_COLS = [5, 6, 7]              # 加速度计数据列索引 [x, y, z]
    
    # 数据格式和单位配置
    TIME_FORMAT = 'ms'                  # 时间格式：'s'(秒) 或 'ms'(毫秒)
    SAMPLE_INTERVAL = 0.01              # IMU采样间隔（秒）
    DATA_FORMAT = 'rate'                # 数据格式：'rate'(速率)或'increment'(增量)
    GYRO_UNIT = 'deg/s'                 # 陀螺仪输入单位：'deg', 'rad', 'deg/s', 'rad/s'
    ACCEL_UNIT = 'g'                    # 加速度计输入单位：'m/s^2' 或 'g'
    
    # 文件读取配置
    SKIP_ROWS = 4300000                       # 跳过的行数（表头等）
    END_ROWS = 5380000                         # 结束行数（可选，留空读取到文件末尾）
    
    # 图片保存配置
    SAVE_IMAGES = False                 # 设置为True保存图像，False不保存
    OUTPUT_FILE = ''                    # 自定义输出文件路径，留空自动生成
    # ================= 配置结束 =================
    
    # 检查输入文件是否存在
    if not os.path.exists(INPUT_FILE):
        print(f"错误: 输入文件 '{INPUT_FILE}' 不存在")
        return 1
    
    # 检查列索引配置是否有效
    if TIME_COL < 0:
        print(f"错误: 时间列索引必须为非负数")
        return 1
    
    if len(GYRO_COLS) != 3:
        print(f"错误: 陀螺仪数据列索引必须为3个值")
        return 1
    
    if len(ACCEL_COLS) != 3:
        print(f"错误: 加速度计数据列索引必须为3个值")
        return 1
    
    # 检查数据格式是否有效
    if DATA_FORMAT not in ['rate', 'increment']:
        print(f"错误: 无效的数据格式 '{DATA_FORMAT}'")
        print("可用选项: 'rate'(速率) 或 'increment'(增量)")
        return 1
    
    # 检查单位是否有效
    valid_gyro_units = ['deg', 'rad', 'deg/s', 'rad/s']
    if GYRO_UNIT not in valid_gyro_units:
        print(f"错误: 无效的陀螺仪单位 '{GYRO_UNIT}'")
        print(f"可用选项: {valid_gyro_units}")
        return 1
    
    valid_accel_units = ['m/s^2', 'g']
    if ACCEL_UNIT not in valid_accel_units:
        print(f"错误: 无效的加速度计单位 '{ACCEL_UNIT}'")
        print(f"可用选项: {valid_accel_units}")
        return 1
    
    # 读取数据
    print("正在读取IMU数据...")
    time_data, gyro_data, accel_data = read_imu_data(
        INPUT_FILE, TIME_COL, GYRO_COLS, 
        ACCEL_COLS, TIME_FORMAT, 
        SKIP_ROWS, END_ROWS-1
    )
    
    if time_data is None:
        print("数据读取失败")
        return 1
    
    print(f"数据读取成功: {len(time_data)} 个数据点")
    print(f"时间数据范围: {time_data.min():.3f} - {time_data.max():.3f} 秒")
    
    # 转换陀螺仪数据
    gyro_converted = convert_gyro_units2rad(
        gyro_data, GYRO_UNIT, DATA_FORMAT, SAMPLE_INTERVAL
    )
    
    # 转换加速度计数据
    accel_converted = convert_accel_units2mg(accel_data, ACCEL_UNIT)
    
    # 构建IMU数据矩阵
    imu_data = np.column_stack([
        gyro_converted,  # 陀螺仪数据 (rad)
        accel_converted,  # 加速度计数据 (mg)
        time_data  # 时间戳 (s)
    ])
    
    # 调用avarimu函数计算Allan方差
    print("开始计算Allan方差...")
    try:
        sigma, tau = avarimu(imu_data, SAMPLE_INTERVAL)
        print("Allan方差计算完成")
        
        # 生成输出文件名（如果未指定）
        if OUTPUT_FILE:
            output_path = OUTPUT_FILE
        else:
            input_stem = Path(INPUT_FILE).stem
            output_path = f"{input_stem}_allan_variance.png"
        
        # 保存图像
        if SAVE_IMAGES:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            print(f"✓ Allan方差图已保存至: {output_path}")
        
        # 显示图像
        plt.show()
        
        print("✓ Allan方差分析完成!")
        
    except Exception as e:
        print(f"Allan方差计算错误: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
        main()