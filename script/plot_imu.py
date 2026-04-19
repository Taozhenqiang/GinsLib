"""
IMU原始数据可视化脚本
功能：读取IMU原始数据，进行单位转换和格式转换，生成六轴数据可视化图表
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys
from pathlib import Path

def main():
    """主函数"""
    # ================= 配置区域 =================
    # 在这里直接设置IMU数据文件的参数配置
    
    # 文件路径配置
    INPUT_FILE = './GNSS/LG69T_Vehicle_complex_20250414/ASM330.txt'  # 输入文件路径
    
    # 数据列索引配置
    TIME_COL = 1                        # 时间列索引（从0开始）
    GYRO_COLS = [2, 3, 4]               # 陀螺仪数据列索引 [x, y, z]
    ACCEL_COLS = [5, 6, 7]              # 加速度计数据列索引 [x, y, z]
    
    # 数据格式和单位配置
    SAMPLE_INTERVAL = 0.01              # IMU采样间隔（秒）
    DATA_FORMAT = 'rate'                # 数据格式：'rate'(速率) 或 'increment'(增量)
    GYRO_UNIT = 'deg/s'                 # 陀螺仪输入单位：'deg', 'rad', 'deg/s', 'rad/s'
    ACCEL_UNIT = 'm/s^2'                # 加速度计输入单位：'m/s^2' 或 'g'
    
    # 文件读取配置
    SKIP_ROWS = 3000                    # 跳过的行数（表头等）
    
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
        ACCEL_COLS, SKIP_ROWS
    )
    
    if time_data is None:
        print("数据读取失败")
        return 1
    
    print(f"数据读取成功: {len(time_data)} 个数据点")
    print(f"时间数据范围: {time_data.min():.3f} - {time_data.max():.3f} 秒")
    print(f"陀螺仪数据形状: {gyro_data.shape}")
    print(f"加速度计数据形状: {accel_data.shape}")
    
    # 数据转换
    print("正在进行数据转换...")
    
    # 转换陀螺仪数据
    gyro_converted = convert_gyro_units(
        gyro_data, GYRO_UNIT, DATA_FORMAT, SAMPLE_INTERVAL
    )
    
    # 转换加速度计数据
    accel_converted = convert_accel_units(accel_data, ACCEL_UNIT)
    
    print("数据转换完成")
    
    # 生成输出文件名（如果未指定）
    if OUTPUT_FILE:
        output_path = OUTPUT_FILE
    else:
        input_stem = Path(INPUT_FILE).stem
        output_path = f"{input_stem}_imu_plot.png"
    
    # 绘制图表
    plot_imu_data(time_data, gyro_converted, accel_converted, 
                  output_path if SAVE_IMAGES else None)
    
    if SAVE_IMAGES:
        print(f"✓ IMU数据可视化完成! 图表已保存至: {output_path}")
    else:
        print("✓ IMU数据可视化完成!")
    
    return 0

def read_imu_data(file_path, time_col, gyro_cols, accel_cols, skip_rows=0):
    """
    读取IMU数据文件
    args:
        file_path: 输入文件路径
        time_col: 时间列索引
        gyro_cols: 陀螺仪数据列索引列表 [x, y, z]
        accel_cols: 加速度计数据列索引列表 [x, y, z]
        skip_rows: 跳过的行数
    return: 时间戳数组，陀螺仪数据数组，加速度计数据数组
    """
    try:
        # 根据文件扩展名选择读取方法
        file_ext = Path(file_path).suffix.lower()
        
        if file_ext in ['.csv']:
            df = pd.read_csv(file_path, skiprows=skip_rows, header=None)
        elif file_ext in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path, skiprows=skip_rows, header=None)
        else:
            # 文本文件处理
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # 跳过指定行数
            data_lines = lines[skip_rows:]
            
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
        
        # 提取时间数据（只提取秒数）
        time_data = df.iloc[:, time_col].values
        
        # 提取陀螺仪数据
        gyro_data = df.iloc[:, gyro_cols].values
        
        # 提取加速度计数据
        accel_data = df.iloc[:, accel_cols].values
        
        return time_data, gyro_data, accel_data
        
    except Exception as e:
        print(f"读取文件错误: {e}")
        return None, None, None

def convert_gyro_units(gyro_data, input_unit, data_format, sample_interval):
    """
    转换陀螺仪数据单位
    args:
        gyro_data: 陀螺仪原始数据
        input_unit: 输入单位 ('deg', 'rad', 'deg/s', 'rad/s')
        data_format: 数据格式 ('rate' 或 'increment')
        sample_interval: 采样间隔（秒）
    return: 转换后的陀螺仪数据（单位：deg/s）
    """
    converted_data = gyro_data.copy()
    
    # 单位转换
    if input_unit == 'deg':
        # 度转换为度/秒（如果已经是速率形式，则不需要转换）
        if data_format == 'increment':
            converted_data = converted_data / sample_interval  
    
    elif input_unit == 'rad':
        # 弧度转换为度/秒
        converted_data = np.degrees(converted_data)
        if data_format == 'increment':
            converted_data = converted_data / sample_interval
    
    elif input_unit == 'rad/s':
        # 弧度/秒转换为度/秒
        converted_data = np.degrees(converted_data)
    
    elif input_unit == 'deg/s':
        # 已经是目标单位，不需要转换
        pass
    
    else:
        print(f"警告: 未知的陀螺仪单位 '{input_unit}'，假设为 deg/s")
    
    return converted_data

def convert_accel_units(accel_data, input_unit):
    """
    转换加速度计数据单位
    args:
        accel_data: 加速度计原始数据
        input_unit: 输入单位 ('m/s^2' 或 'g')
    return: 转换后的加速度计数据（单位：g）
    """
    converted_data = accel_data.copy()
    
    if input_unit == 'm/s^2':
        # m/s^2 转换为 g (1 g = 9.80665 m/s^2)
        converted_data = converted_data / 9.80665
    
    elif input_unit == 'g':
        # 已经是目标单位，不需要转换
        pass
    
    else:
        print(f"警告: 未知的加速度计单位 '{input_unit}'，假设为 g")
    
    return converted_data

def plot_imu_data(time_data, gyro_data, accel_data, output_file=None):
    """
    绘制IMU六轴数据
    args:
        time_data: 时间戳数据
        gyro_data: 陀螺仪数据（deg/s）
        accel_data: 加速度计数据（g）
        output_file: 输出文件路径（可选）
    """
    # 创建图形和子图
    fig, axes = plt.subplots(3, 2, figsize=(15, 10))
    
    # 陀螺仪数据标签
    gyro_labels = ['Gyro X (deg/s)', 'Gyro Y (deg/s)', 'Gyro Z (deg/s)']
    gyro_colors = ['red', 'green', 'blue']
    
    # 加速度计数据标签
    accel_labels = ['Accel X (g)', 'Accel Y (g)', 'Accel Z (g)']
    accel_colors = ['red', 'green', 'blue']
    
    # 绘制陀螺仪数据（左侧三个子图）
    for i in range(3):
        ax = axes[i, 0]
        ax.plot(time_data, gyro_data[:, i], color=gyro_colors[i], linewidth=1)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel(gyro_labels[i])
        ax.grid(True, alpha=0.3)
    
    # 绘制加速度计数据（右侧三个子图）
    for i in range(3):
        ax = axes[i, 1]
        ax.plot(time_data, accel_data[:, i], color=accel_colors[i], linewidth=1)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel(accel_labels[i])
        ax.grid(True, alpha=0.3)
    
    # 调整布局
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    
    # 保存或显示图形
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"图表已保存至: {output_file}")
    
    plt.show()


if __name__ == "__main__":
    main()