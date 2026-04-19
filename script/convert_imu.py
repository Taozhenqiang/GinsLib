import numpy as np
import os
import re
import pandas as pd
import argparse
import sys
from time_conversion import epoch2time, utc2gpst, time2gpst, GTimeT
from common import cal_gravity

# ========================= 配置区域 =========================
# 在这里修改配置参数（如果使用命令行参数，这些配置将被覆盖）

# 输入文件路径
INPUT_FILE = "./GNSS/LG69T_Vehicle_complex_20250414/LG69TAP01-TEMP0711-1HZ_INS.dat"  

# 要跳过的行数（通常是表头或注释行）
SKIP_LINES = 1 # 5000

# 时间系统配置
TIME_SYSTEM = "UNIX"  # 可选: "GPST", "UTC", "UNIX"
INPUT_TIME_FORMAT = "unix"  # 可选: "calendar"(年月日时分秒, 6列), "week_tow"(周与周内秒, 2列), "unix"(UNIX时间, 1列)
OUTPUT_TIME_FORMAT = "week_tow"  # 输出格式固定为周/周内秒

# IMU数据列索引（从0开始计数）
GRO_COLS = [3, 4, 5]
ACC_COLS = [6, 7, 8]
IMU_COLS = GRO_COLS + ACC_COLS  # [gyro_x, gyro_y, gyro_z, acc_x, acc_y, acc_z]

ACC_UNIT = "g" # 加速度计的单位是 m/s^2 还是 g
LAT_HEIGHT = [31.81845, 44] # 所在地区的纬度与高程

# 时间戳列索引（根据INPUT_TIME_FORMAT确定）
TIME_COLS = [1]  # 默认: 年月日时分秒格式 # [0, 1, 2, 3, 4, 5]

# IMU采样间隔（秒）
IMU_SAMPLE_INTERVAL = 0.01  # 默认10ms采样间隔

# 输出文件路径（如果为None，则自动生成）
OUTPUT_FILE = './GNSS/LG69T_Vehicle_complex_20250414/asm330lhh.txt'

# 日志文件路径
LOG_FILE = './GNSS/LG69T_Vehicle_complex_20250414/asm330lhh.log'  # None表示自动生成

# 文件分隔符（手动设置，避免自动检测错误）
DELIMITER = ''  # 建议手动指定，如",", "\t", " "等

# 目标指令过滤（只保留包含这些指令的行）
TARGET_COMMANDS = ["$PQTMRAWIMU"]  # 示例：只保留包含$PQTMRAWIMU的行
# TARGET_COMMANDS = ["$PQTMRAWIMU", "$GPGGA"]  # 示例：保留多种指令
# TARGET_COMMANDS = []  # 空列表表示不过滤，保留所有行

# 是否显示详细信息
VERBOSE = True

# 输出格式配置
FIELD_WIDTH = 15  # 每个字段的宽度（字符数）
DECIMAL_PLACES = 6  # 小数位数

# ========================= 函数定义 =========================

def detect_delimiter(line):
    """
    检测行中的分隔符
    """
    # 统计各种分隔符的数量
    space_count = len(re.findall(r'\s+', line))
    comma_count = line.count(',')
    tab_count = line.count('\t')
    semicolon_count = line.count(';')
    asterisk_count = line.count('*')
    
    # 返回出现次数最多的分隔符
    counts = {
        'space': space_count,
        'comma': comma_count, 
        'tab': tab_count,
        'semicolon': semicolon_count,
        'asterisk': asterisk_count
    }
    
    if max(counts.values()) > 0:
        return max(counts, key=counts.get)
    else:
        return 'space'  # 默认使用空格

def parse_line(line):
    """
    根据分隔符解析行数据，支持多种分隔符组合
    """
    # 使用正则表达式处理多种分隔符
    pattern = r'[\s,;\*]+'  # 空格、逗号、分号、星号
    
    # 使用正则表达式分割
    elements = re.split(pattern, line.strip())
    
    # 过滤空字符串
    return [elem for elem in elements if elem]

def filter_target_lines(lines, target_commands):
    """
    过滤行，只保留包含目标指令的行
    """
    if not target_commands:  # 空列表表示不过滤
        return lines
    
    filtered_lines = []
    for line in lines:
        # 检查行是否包含任何目标指令
        if any(cmd in line for cmd in target_commands):
            filtered_lines.append(line)
    
    return filtered_lines

def unix_time_to_gpst(unix_time):
    """
    将UNIX时间转换为GPS时间
    UNIX时间: 1970年1月1日以来的秒数
    参考time_conversion.py的实现方式
    """
    # UNIX时间起点: 1970-01-01 00:00:00
    
    # 计算GPS时间
    t_gps = GTimeT()
    t_gps.time = int(unix_time)
    t_gps.sec = unix_time - int(unix_time)
    
    # 转换为周和周内秒
    week = [0]
    tow = time2gpst(t_gps, week)
    
    return week[0], tow

def convert_time_format(time_data, input_format, time_system):
    """
    转换时间格式为周/周内秒
    """
    converted_time = []
    
    for i, time_row in enumerate(time_data):
        try:
            if input_format == "calendar":
                # 输入为格里高利历格式 [年, 月, 日, 时, 分, 秒]
                if len(time_row) >= 6:
                    # 检查数据有效性
                    if not all(isinstance(x, (int, float)) for x in time_row[:6]):
                        if VERBOSE:
                            print(f"警告: 第{i+1}行包含非数字时间数据: {time_row[:6]}")
                        converted_time.append([0.0, 0.0])
                        continue
                    
                    # 将日历时间转换为计算机时
                    t_epoch = epoch2time(time_row[:6])
                    
                    # 检查时间转换是否有效
                    if t_epoch.time == 0:
                        if VERBOSE:
                            print(f"警告: 第{i+1}行时间转换失败: {time_row[:6]}")
                        converted_time.append([0.0, 0.0])
                        continue
                    
                    # 根据时间系统转换为GPS时间
                    if time_system == "UTC":
                        t_gps = utc2gpst(t_epoch)
                    else:  # GPST
                        t_gps = t_epoch
                    
                    # 转换为周和周内秒
                    week = [0]
                    tow = time2gpst(t_gps, week)
                    
                    # 检查转换结果
                    if week[0] == 0 and tow == 0:
                        if VERBOSE:
                            print(f"警告: 第{i+1}行GPS时间转换异常: {time_row[:6]} -> 周:{week[0]}, 秒:{tow}")
                        converted_time.append([0.0, 0.0])
                    else:
                        converted_time.append([week[0], tow])
                else:
                    if VERBOSE:
                        print(f"警告: 第{i+1}行数据不足6个时间元素: {time_row}")
                    converted_time.append([0.0, 0.0])
                    
            elif input_format == "week_tow":
                # 输入已经是周/周内秒格式
                if len(time_row) >= 2:
                    # 检查数据有效性
                    if not all(isinstance(x, (int, float)) for x in time_row[:2]):
                        if VERBOSE:
                            print(f"警告: 第{i+1}行包含非数字周/秒数据: {time_row[:2]}")
                        converted_time.append([0.0, 0.0])
                        continue
                    
                    converted_time.append(time_row[:2])
                else:
                    if VERBOSE:
                        print(f"警告: 第{i+1}行数据不足2个时间元素: {time_row}")
                    converted_time.append([0.0, 0.0])
                    
            elif input_format == "unix":
                # 输入为UNIX时间格式
                if len(time_row) >= 1:
                    # 检查数据有效性
                    if not isinstance(time_row[0], (int, float)):
                        if VERBOSE:
                            print(f"警告: 第{i+1}行包含非数字UNIX时间: {time_row[0]}")
                        converted_time.append([0.0, 0.0])
                        continue
                    
                    # 转换为GPS周和秒
                    week, tow = unix_time_to_gpst(time_row[0])
                    converted_time.append([week, tow])
                else:
                    if VERBOSE:
                        print(f"警告: 第{i+1}行数据不足1个时间元素: {time_row}")
                    converted_time.append([0.0, 0.0])
                    
            else:
                if VERBOSE:
                    print(f"错误: 第{i+1}行未知输入时间格式: {input_format}")
                converted_time.append([0.0, 0.0])
                
        except Exception as e:
            if VERBOSE:
                print(f"错误: 第{i+1}行时间转换异常 - 行内容: {time_row}, 错误: {e}")
            converted_time.append([0.0, 0.0])
    
    # 统计转换结果
    valid_count = sum(1 for item in converted_time if item != [0.0, 0.0])
    invalid_count = len(converted_time) - valid_count
    
    if VERBOSE:
        if invalid_count > 0:
            print(f"时间格式转换完成: 有效{valid_count}行, 无效{invalid_count}行")
        else:
            print(f"时间格式转换完成: 全部{valid_count}行有效")
    
    return np.array(converted_time)

def read_csv_file(file_path, skip_lines=0, delimiter=None, target_commands=None):
    """读取CSV文件并返回所有数据"""
    try:
        # 先读取所有行进行过滤
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # 过滤目标指令行
        if target_commands:
            original_count = len(lines)
            lines = filter_target_lines(lines, target_commands)
            filtered_count = len(lines)
            if VERBOSE:
                print(f"指令过滤: 原始{original_count}行 -> 保留{filtered_count}行")
        
        # 跳过指定行数
        if skip_lines >= len(lines):
            if VERBOSE:
                print(f"警告: 跳过的行数({skip_lines})超过文件总行数({len(lines)})")
            skip_lines = 0
        
        lines = lines[skip_lines:]
        
        if delimiter:
            # 使用指定分隔符
            # 将过滤后的行写入临时字符串，然后读取
            from io import StringIO
            content = ''.join(lines)
            df = pd.read_csv(StringIO(content), header=None, sep=delimiter, dtype=str)
        else:
            # 自动检测分隔符
            if lines:
                detected_delimiter = detect_delimiter(lines[0])
                if detected_delimiter == 'space':
                    sep = r'\s+'
                elif detected_delimiter == 'comma':
                    sep = ','
                elif detected_delimiter == 'tab':
                    sep = '\t'
                elif detected_delimiter == 'semicolon':
                    sep = ';'
                else:
                    sep = None
            else:
                sep = None
            
            if sep:
                from io import StringIO
                content = ''.join(lines)
                df = pd.read_csv(StringIO(content), header=None, sep=sep, dtype=str, engine='python')
            else:
                df = pd.DataFrame()
        
        if df.empty:
            return np.array([])
        
        # 将所有数据转换为浮点数，无法转换的设为NaN
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # 检查是否有NaN值
        nan_count = df.isna().sum().sum()
        if nan_count > 0 and VERBOSE:
            print(f"警告: CSV文件中发现{nan_count}个无法转换为数字的值，已设为NaN")
        
        # 将NaN替换为0
        df = df.fillna(0)
        
        return df.to_numpy()
    except Exception as e:
        print(f"读取CSV文件错误: {e}")
        return None

def read_excel_file(file_path, skip_lines=0, target_commands=None):
    """读取Excel文件并返回所有数据"""
    try:
        df = pd.read_excel(file_path, skiprows=skip_lines, header=None, dtype=str)
        
        # Excel文件不支持行级过滤，但可以检查第一列是否包含目标指令
        if target_commands and not df.empty:
            original_count = len(df)
            # 假设指令在第一列
            mask = df.iloc[:, 0].apply(lambda x: any(cmd in str(x) for cmd in target_commands))
            df = df[mask]
            filtered_count = len(df)
            if VERBOSE:
                print(f"指令过滤: 原始{original_count}行 -> 保留{filtered_count}行")
        
        # 将所有数据转换为浮点数，无法转换的设为NaN
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # 检查是否有NaN值
        nan_count = df.isna().sum().sum()
        if nan_count > 0 and VERBOSE:
            print(f"警告: Excel文件中发现{nan_count}个无法转换为数字的值，已设为NaN")
        
        # 将NaN替换为0
        df = df.fillna(0)
        
        return df.to_numpy()
    except Exception as e:
        print(f"读取Excel文件错误: {e}")
        return None

def read_text_file(file_path, skip_lines=0, target_commands=None):
    """读取文本文件并返回所有数据"""
    try:
        # 逐行读取，跳过编码错误的行
        lines = []
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f_in:
            for line_num, line in enumerate(f_in, 1):
                lines.append(line)
        
        if VERBOSE:
            print(f"成功读取 {len(lines)} 行数据")
        
        # 过滤目标指令行
        if target_commands:
            original_count = len(lines)
            lines = filter_target_lines(lines, target_commands)
            filtered_count = len(lines)
            if VERBOSE:
                print(f"指令过滤: 原始{original_count}行 -> 保留{filtered_count}行")
        
        # 跳过指定行数
        if skip_lines >= len(lines):
            if VERBOSE:
                print(f"警告: 跳过的行数({skip_lines})超过文件总行数({len(lines)})")
            skip_lines = 0
        
        data_lines = lines[skip_lines:]
        
        
        # 处理数据
        processed_data = []
        
        for i, line in enumerate(data_lines, skip_lines + 1):
            line = line.strip()
            
            # 跳过空行和注释行
            if not line or line.startswith('#') or line.startswith('%'):
                continue
            
            # 解析行数据
            try:
                elements = parse_line(line)
                
                # 将非数字元素替换为占位符，保持每行长度一致
                numeric_elements = []
                for elem in elements:
                    try:
                        numeric_elements.append(float(elem))
                    except ValueError:
                        # 非数字元素替换为占位符999
                        numeric_elements.append(999.0)
                
                processed_data.append(numeric_elements)
                
            except Exception as e:
                if VERBOSE:
                    print(f"警告: 第{i}行解析失败: {e}")
                continue
        
        # 检查数据长度一致性并转换为numpy数组
        if not processed_data:
            return None
        
        # 统计每行的长度
        row_lengths = [len(row) for row in processed_data]
        
        # 找到最常见的行长度（众数）
        from collections import Counter
        length_counter = Counter(row_lengths)
        most_common_length, most_common_count = length_counter.most_common(1)[0]
        
        if VERBOSE:
            print(f"数据行长度统计: {dict(length_counter)}")
        
        # 剔除长度异常的行
        valid_data = []
        abnormal_indices = []
        
        for i, row in enumerate(processed_data):
            if len(row) == most_common_length:
                valid_data.append(row)
            else:
                abnormal_indices.append(row[1])
        
        if abnormal_indices and VERBOSE:
            print(f"剔除 {len(abnormal_indices)} 行长度异常的数据")
            if len(abnormal_indices) <= 10:  # 只显示前10个异常行
                print(f"异常行时间索引: {abnormal_indices[:10]}")
            else:
                print(f"异常行时间索引: {abnormal_indices[:10]} ... (共{len(abnormal_indices)}行)")
        
        if not valid_data:
            print("错误: 没有找到有效的数据行")
            return None
        
        # 转换为numpy数组
        data_array = np.array(valid_data)
        
        return data_array
    except Exception as e:
        print(f"读取文本文件错误: {e}")
        return None

def check_imu_data_quality(time_data, imu_data, sample_interval, log_file=None):
    """
    检查IMU数据质量，统计时间间隔异常
    """
    if len(time_data) < 2:
        return [], 0, 0
    
    abnormal_indices = []
    missing_count = 0
    
    # 计算时间间隔
    time_intervals = []
    for i in range(1, len(time_data)):
        # 计算相邻时间戳的差值（秒）
        time_diff = (time_data[i][0] - time_data[i-1][0]) * 7 * 24 * 3600 + (time_data[i][1] - time_data[i-1][1])
        time_intervals.append(time_diff)
        
        # 检查时间间隔是否异常
        if abs(time_diff - sample_interval) > sample_interval * 0.1:  # 允许10%的误差
            abnormal_indices.append(i)
            if time_diff > sample_interval * 1.5:  # 大于1.5倍采样间隔视为缺失
                missing_count += int(time_diff / sample_interval) - 1
    
    # 写入日志文件
    if log_file:
        with open(log_file, 'w', encoding='utf-8') as f:
            f.write("IMU数据质量检查报告\n")
            f.write("=" * 50 + "\n")
            f.write(f"总数据行数: {len(time_data)}\n")
            f.write(f"采样间隔: {sample_interval}秒\n")
            f.write(f"异常数据行数: {len(abnormal_indices)}\n")
            f.write(f"估计缺失数据条数: {missing_count}\n")
            f.write("\n异常数据详情:\n")
            
            for idx in abnormal_indices:
                time_diff = time_intervals[idx-1]
                f.write(f"第{idx+1}行: 时间间隔={time_diff:.6f}秒 (期望={sample_interval}秒)\n")
            
            f.write("\n时间间隔统计:\n")
            if time_intervals:
                f.write(f"最小时间间隔: {min(time_intervals):.6f}秒\n")
                f.write(f"最大时间间隔: {max(time_intervals):.6f}秒\n")
                f.write(f"平均时间间隔: {np.mean(time_intervals):.6f}秒\n")
                f.write(f"时间间隔标准差: {np.std(time_intervals):.6f}秒\n")
    
    return abnormal_indices, missing_count, time_intervals

def generate_right_aligned_format(num_columns, field_width=15, decimal_places=6):
    """
    生成右对齐的格式化字符串
    """
    format_strings = [f"%{field_width}.{decimal_places}f" for _ in range(num_columns)]
    return ' '.join(format_strings)

def extract_imu_data():
    """
    提取IMU数据并转换为标准格式
    """
    # 使用全局配置
    input_file = INPUT_FILE
    skip_lines = SKIP_LINES
    time_cols = TIME_COLS
    imu_cols = IMU_COLS
    acc_unit = ACC_UNIT
    time_system = TIME_SYSTEM
    input_time_format = INPUT_TIME_FORMAT
    output_file = OUTPUT_FILE
    log_file = LOG_FILE
    sample_interval = IMU_SAMPLE_INTERVAL
    delimiter = DELIMITER
    target_commands = TARGET_COMMANDS
    verbose = VERBOSE
    field_width = FIELD_WIDTH
    decimal_places = DECIMAL_PLACES
    
    # 设置默认输出文件名
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_imu.truth"
    
    # 设置默认日志文件名
    if log_file is None:
        base_name = os.path.splitext(input_file)[0]
        log_file = f"{base_name}_imu_quality.log"
    
    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在")
        return False
    
    try:
        # 根据文件扩展名选择读取方法
        file_ext = os.path.splitext(input_file)[1].lower()
        
        if file_ext == '.csv':
            data_array = read_csv_file(input_file, skip_lines, delimiter, target_commands)
        elif file_ext in ['.xlsx', '.xls']:
            data_array = read_excel_file(input_file, skip_lines, target_commands)
        else:
            data_array = read_text_file(input_file, skip_lines, target_commands)
        
        if data_array is None:
            print("错误: 无法读取文件数据")
            return False
        
        # 提取时间数据
        time_data = []
        for row in data_array:
            if len(row) > max(time_cols):
                time_row = [row[col] for col in time_cols if col < len(row)]
                time_data.append(time_row)
            else:
                time_data.append([])
        
        # 转换时间格式
        converted_time = convert_time_format(time_data, input_time_format, time_system)
        
        # 提取IMU数据
        imu_data = []
        for row in data_array:
            if len(row) > max(imu_cols):
                imu_row = [row[col] for col in imu_cols if col < len(row)]
                # 转换加速度计单位为 m/s^2
                if acc_unit == 'g':
                    g = cal_gravity(np.deg2rad(LAT_HEIGHT[0]), LAT_HEIGHT[1])
                    imu_row[3] = imu_row[3] * g
                    imu_row[4] = imu_row[4] * g
                    imu_row[5] = imu_row[5] * g
                imu_data.append(imu_row)
            else:
                imu_data.append([0.0] * len(imu_cols))
        
        imu_array = np.array(imu_data)
        
        # 合并时间和IMU数据
        if len(converted_time) != len(imu_array):
            print(f"错误: 时间数据({len(converted_time)}行)和IMU数据({len(imu_array)}行)行数不匹配")
            return False
        
        # 只保留有效的时间数据
        valid_indices = [i for i, time_row in enumerate(converted_time) if not np.array_equal(time_row, [0.0, 0.0])]
        
        if not valid_indices:
            print("错误: 没有找到有效的时间数据")
            return False
        
        # 提取有效的时间数据用于质量检查
        valid_time_data = [converted_time[i] for i in valid_indices]
        valid_imu_data = [imu_array[i] for i in valid_indices]
        
        # 检查IMU数据质量
        abnormal_indices, missing_count, time_intervals = check_imu_data_quality(
            valid_time_data, valid_imu_data, sample_interval, log_file)
        
        if verbose:
            print(f"IMU数据质量检查完成:")
            print(f"  有效数据行数: {len(valid_time_data)}")
            print(f"  异常数据行数: {len(abnormal_indices)}")
            print(f"  估计缺失数据条数: {missing_count}")
            if time_intervals:
                print(f"  时间间隔统计: 最小={min(time_intervals):.6f}s, 最大={max(time_intervals):.6f}s, 平均={np.mean(time_intervals):.6f}s")
        
        # 构建最终数据（包含所有有效数据）
        final_data = []
        for i in range(len(valid_time_data)):
            row = np.concatenate([valid_time_data[i], valid_imu_data[i]])
            final_data.append(row)
        
        final_array = np.array(final_data)
        
        if verbose:
            print(f"最终数据维度: {final_array.shape}")
        
        # 生成右对齐的格式化字符串
        num_columns = final_array.shape[1]
        fmt_string = generate_right_aligned_format(num_columns, field_width, decimal_places)
        
        # 保存为标准IMU文件（右对齐格式）
        # 注意：覆盖源文件时需要先关闭文件，否则会报错
        np.savetxt(output_file, final_array, fmt=fmt_string, delimiter=' ')
        
        print(f"\n转换完成!")
        
        # 显示前几行数据预览
        if verbose and len(final_array) > 0:
            print(f"\n数据预览（前3行，右对齐格式）:")
            for i in range(min(3, len(final_array))):
                formatted_line = fmt_string % tuple(final_array[i])
                print(f"第{i+1}行: {formatted_line}")
        
        return True
        
    except Exception as e:
        print(f"转换过程中发生错误: {e}")
        return False

def main():
    """
    主函数
    """
    
    print("开始提取IMU数据...")
    print("=" * 50)
    
    # 显示当前配置
    print(f"输入文件: {INPUT_FILE}")
    print(f"输出文件: {OUTPUT_FILE}")
    print(f"日志文件: {LOG_FILE}")
    print(f"跳过行数: {SKIP_LINES}")
    print(f"采样间隔: {IMU_SAMPLE_INTERVAL}秒")
    print(f"时间系统: {TIME_SYSTEM}")
    print(f"输入时间格式: {INPUT_TIME_FORMAT}")
    print(f"时间列索引: {TIME_COLS}")
    print(f"IMU列索引: {IMU_COLS}")
    print(f"采样间隔: {IMU_SAMPLE_INTERVAL}秒")
    print(f"目标指令: {TARGET_COMMANDS}")
    print("=" * 50)
    
    # 执行转换
    success = extract_imu_data()
    
    if success:
        print("\nIMU数据提取成功完成!")
    else:
        print("\nIMU数据提取失败!")

if __name__ == "__main__":
    main()