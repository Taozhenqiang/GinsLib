import numpy as np
import os
import re
import pandas as pd
import csv
from coordinate_transformation import pos2ecef
from coordinate_transformation import xyz2enu
from time_conversion import epoch2time, utc2gpst, time2gpst, GTimeT

# ========================= 配置区域 =========================
# 在这里修改配置参数

# 输入文件路径
INPUT_FILE = "./GNSS/LG69T_Vehicle_complex_20250414/POS.tru"  

# 要跳过的行数（通常是表头或注释行）
SKIP_LINES = 108

# 时间系统配置
TIME_SYSTEM = "GPST"  # 可选: "GPST" 或 "UTC"
INPUT_TIME_FORMAT = "week_tow"  # 可选: "week_tow" 或 "calendar"
OUTPUT_TIME_FORMAT = "week_tow"  # 输出格式固定为周/周内秒

# 单独输入各数据项的列索引（从0开始计数，如果某项不存在则设为-1）
# 时间戳列索引（必须提供）
TIME_COLS = [0, 1]  # [周, 秒] 或 [年, 月, 日, 时, 分, 秒] 根据INPUT_TIME_FORMAT确定

# 位置相关列索引（如果存在）
POSITION_COLS = [2, 3, 4]  # [x, y, z] 或 [lat, lon, height]
POSITION_FORMAT = "geodetic"  # 可选: "ecef" 或 "geodetic"

# 速度相关列索引（如果存在）
VELOCITY_COLS = [6, 5, 7]  # [vx, vy, vz] 或 [ve, vn, vu]
VELOCITY_FORMAT = "enu"  # 可选: "ecef" 或 "enu"
VELOCITY_COEFFICIENTS = [1, 1, -1]  # 速度系数 [x系数, y系数, z系数]

# 姿态相关列索引（如果存在）
ATTITUDE_COLS = [9, 8, 10]  # [pitch, roll, heading]

# 输出文件路径（如果为None，则自动生成）
OUTPUT_FILE = "./GNSS/LG69T_Vehicle_complex_20250414/truth_pva.truth"

# 是否显示详细信息
VERBOSE = True

# 输出格式配置
FIELD_WIDTH = 15  # 每个字段的宽度（字符数）
DECIMAL_PLACES = 6  # 小数位数

# ========================= 函数定义 =========================

def detect_delimiter(line):
    """检测行中的分隔符"""
    space_count = len(re.findall(r'\s+', line))
    comma_count = line.count(',')
    tab_count = line.count('\t')
    
    if space_count > comma_count and space_count > tab_count:
        return 'space'
    elif comma_count > space_count and comma_count > tab_count:
        return 'comma'
    elif tab_count > space_count and tab_count > comma_count:
        return 'tab'
    else:
        return 'space'

def parse_line(line, delimiter):
    """根据分隔符解析行数据"""
    if delimiter == 'space':
        return re.split(r'\s+', line.strip())
    elif delimiter == 'comma':
        return line.strip().split(',')
    elif delimiter == 'tab':
        return line.strip().split('\t')
    else:
        return line.strip().split()

def generate_right_aligned_format(num_columns, field_width=15, decimal_places=6):
    """生成右对齐的格式化字符串"""
    format_strings = [f"%{field_width}.{decimal_places}f" for _ in range(num_columns)]
    return ' '.join(format_strings)

def convert_geodetic_to_ecef(position_data):
    """将大地坐标转换为ECEF坐标"""
    ecef_positions = []
    for pos in position_data:
        if len(pos) == 3:
            # 假设输入是[纬度, 经度, 高度]（度）
            lat_rad = np.deg2rad(pos[0])
            lon_rad = np.deg2rad(pos[1])
            height = pos[2]
            llh_pos = [lat_rad, lon_rad, height]
            ecef = pos2ecef(llh_pos)
            ecef_positions.append(ecef)
        else:
            ecef_positions.append([0.0, 0.0, 0.0])  # 无效数据
    return np.array(ecef_positions)

def convert_enu_to_ecef(velocity_data, position_data, coefficients):
    """将ENU速度转换为ECEF速度"""
    ecef_velocities = []
    for i, (vel, pos) in enumerate(zip(velocity_data, position_data)):
        if len(vel) == 3 and len(pos) >= 2:
            # 应用速度系数
            vel_scaled = [vel[0] * coefficients[0], vel[1] * coefficients[1], vel[2] * coefficients[2]]
            
            # 计算变换矩阵（只需要纬度和经度）
            E = xyz2enu(np.array(pos[:2])*np.pi/180.0)
            
            # ENU到ECEF的变换：v_ecef = E^T * v_enu
            v_ecef = np.dot(E.T, vel_scaled)
            ecef_velocities.append(v_ecef)
        else:
            ecef_velocities.append([0.0, 0.0, 0.0])  # 无效数据
    return np.array(ecef_velocities)

def convert_time_format(time_data, input_format, time_system, output_format="week_tow"):
    """
    转换时间格式
    args: time_data - 时间数据列表
           input_format - 输入时间格式 ("week_tow" 或 "calendar")
           time_system - 时间系统 ("GPST" 或 "UTC")
           output_format - 输出时间格式 (固定为"week_tow")
    return: 转换后的时间数据 [周, 周内秒]
    """
    converted_time = []
    
    for i, time_row in enumerate(time_data):
        try:
            if input_format == "calendar":
                # 输入为格里高利历格式 [年, 月, 日, 时, 分, 秒]
                if len(time_row) >= 6:
                    # 检查数据有效性
                    if not all(isinstance(x, (int, float)) for x in time_row[:6]):
                        print(f"警告: 第{i+1}行包含非数字时间数据: {time_row[:6]}")
                        converted_time.append([0.0, 0.0])
                        continue
                    
                    # 将日历时间转换为计算机时
                    t_epoch = epoch2time(time_row[:6])
                    
                    # 检查时间转换是否有效
                    if t_epoch.time == 0:
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
                        print(f"警告: 第{i+1}行GPS时间转换异常: {time_row[:6]} -> 周:{week[0]}, 秒:{tow}")
                        converted_time.append([0.0, 0.0])
                    else:
                        converted_time.append([week[0], tow])
                else:
                    print(f"警告: 第{i+1}行数据不足6个时间元素: {time_row}")
                    converted_time.append([0.0, 0.0])  # 无效数据
                    
            elif input_format == "week_tow":
                # 输入已经是周/周内秒格式
                if len(time_row) >= 2:
                    # 检查数据有效性
                    if not all(isinstance(x, (int, float)) for x in time_row[:2]):
                        print(f"警告: 第{i+1}行包含非数字周/秒数据: {time_row[:2]}")
                        converted_time.append([0.0, 0.0])
                        continue
                    
                    converted_time.append(time_row[:2])
                else:
                    print(f"警告: 第{i+1}行数据不足2个时间元素: {time_row}")
                    converted_time.append([0.0, 0.0])  # 无效数据
                    
            else:
                print(f"错误: 第{i+1}行未知输入时间格式: {input_format}")
                converted_time.append([0.0, 0.0])  # 未知格式
                
        except Exception as e:
            print(f"错误: 第{i+1}行时间转换异常 - 行内容: {time_row}, 错误: {e}")
            converted_time.append([0.0, 0.0])
    
    # 统计转换结果
    valid_count = sum(1 for item in converted_time if item != [0.0, 0.0])
    invalid_count = len(converted_time) - valid_count
    
    if invalid_count > 0:
        print(f"时间格式转换完成: 有效{valid_count}行, 无效{invalid_count}行")
    else:
        print(f"时间格式转换完成: 全部{valid_count}行有效")
    
    return np.array(converted_time)

def read_csv_file(file_path, skip_lines=0):
    """读取CSV文件并返回所有数据"""
    try:
        # 强制将所有列转换为数字类型，无法转换的设为NaN
        df = pd.read_csv(file_path, skiprows=skip_lines, header=None, dtype=str)
        
        # 将所有数据转换为浮点数，无法转换的设为NaN
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # 检查是否有NaN值
        nan_count = df.isna().sum().sum()
        if nan_count > 0:
            print(f"警告: CSV文件中发现{nan_count}个无法转换为数字的值，已设为NaN")
        
        # 将NaN替换为0或删除包含NaN的行
        df = df.fillna(0)  # 或者使用 df = df.dropna() 删除包含NaN的行
        
        return df.to_numpy()
    except Exception as e:
        print(f"读取CSV文件错误: {e}")
        return None

def read_excel_file(file_path, skip_lines=0, sheet_name=0):
    """读取Excel文件并返回所有数据"""
    try:
        # 强制将所有列转换为字符串，然后转换为数字
        df = pd.read_excel(file_path, skiprows=skip_lines, header=None, sheet_name=sheet_name, dtype=str)
        
        # 将所有数据转换为浮点数，无法转换的设为NaN
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # 检查是否有NaN值
        nan_count = df.isna().sum().sum()
        if nan_count > 0:
            print(f"警告: Excel文件中发现{nan_count}个无法转换为数字的值，已设为NaN")
            
            # 显示包含NaN的具体位置
            nan_rows = df.isna().any(axis=1)
            if nan_rows.any():
                print("包含NaN值的行号:", df[nan_rows].index.tolist())
        
        # 将NaN替换为0
        df = df.fillna(0)
        
        return df.to_numpy()
    except Exception as e:
        print(f"读取Excel文件错误: {e}")
        return None

def read_text_file(file_path, skip_lines=0):
    """读取文本文件并返回所有数据"""
    try:
        with open(file_path, 'r', encoding='gbk') as f_in:
            lines = f_in.readlines()
        
        if skip_lines >= len(lines):
            print(f"警告: 跳过的行数({skip_lines})超过文件总行数({len(lines)})")
            skip_lines = 0
        
        data_lines = lines[skip_lines:]
        
        # 检测分隔符
        delimiter = 'space'
        for line in data_lines:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('%'):
                delimiter = detect_delimiter(line)
                break
        
        # 处理数据
        processed_data = []
        
        for line in data_lines:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('%'):
                continue
            
            try:
                elements = parse_line(line, delimiter)
                numeric_elements = []
                for elem in elements:
                    try:
                        numeric_elements.append(float(elem))
                    except ValueError:
                        continue
                
                if numeric_elements:
                    processed_data.append(numeric_elements)
                
            except Exception:
                continue
        
        if not processed_data:
            return None
        
        return np.array(processed_data)
    except Exception as e:
        print(f"读取文本文件错误: {e}")
        return None

def extract_and_convert_data(data_array, time_cols, time_format, time_system, 
                           position_cols, position_format, 
                           velocity_cols, velocity_format, velocity_coefficients, 
                           attitude_cols):
    """提取并转换数据，缺失数据自动置零"""
    if data_array is None or len(data_array) == 0:
        return None
    
    # 计算输出数据的列数
    output_columns = 0
    if time_cols and all(col >= 0 for col in time_cols):
        output_columns +=len(time_cols)  # 时间戳为2列或6列（周,周内秒/格里高利历）
    if position_cols and all(col >= 0 for col in position_cols):
        output_columns += 3  # 位置数据固定为3列（ECEF坐标）
    if velocity_cols and all(col >= 0 for col in velocity_cols):
        output_columns += 3  # 速度数据固定为3列（ECEF速度）
    if attitude_cols and all(col >= 0 for col in attitude_cols):
        output_columns += 3  # 姿态数据固定为3列（pitch, roll, heading）
    
    # 提取时间戳数据
    time_data = []
    has_time_data = time_cols and all(col >= 0 for col in time_cols)
    if has_time_data:
        for row in data_array:
            time_row = []
            for col in time_cols:
                if col < len(row):
                    time_row.append(row[col])
                else:
                    time_row.append(0.0)
            time_data.append(time_row)
        
        # 时间格式转换
        time_data = convert_time_format(time_data, time_format, time_system)
    else:
        # 如果没有时间戳数据，创建全零的时间戳
        for _ in range(len(data_array)):
            time_data.append([0.0, 0.0])  # 默认2列时间戳（周, 周内秒）
    
    # 提取位置数据
    position_data = []
    has_position_data = position_cols and all(col >= 0 for col in position_cols)
    if has_position_data:
        for row in data_array:
            pos_row = []
            for col in position_cols:
                if col < len(row):
                    pos_row.append(row[col])
                else:
                    pos_row.append(0.0)
            position_data.append(pos_row)
        
        # 坐标转换：大地坐标转ECEF
        if position_format.lower() == "geodetic":
            ecef_position = convert_geodetic_to_ecef(position_data)
        else:
            ecef_position = position_data
    else:
        # 如果没有位置数据，创建全零的位置数据
        for _ in range(len(data_array)):
            ecef_position.append([0.0, 0.0, 0.0])
    
    # 提取速度数据
    velocity_data = []
    has_velocity_data = velocity_cols and all(col >= 0 for col in velocity_cols)
    if has_velocity_data:
        for row in data_array:
            vel_row = []
            for col in velocity_cols:
                if col < len(row):
                    vel_row.append(row[col])
                else:
                    vel_row.append(0.0)
            velocity_data.append(vel_row)
        
        # 坐标转换：ENU转ECEF
        if velocity_format.lower() == "enu" and len(position_data) > 0:
            ecef_velocity = convert_enu_to_ecef(velocity_data, position_data, velocity_coefficients)
        else:
            ecef_velocity = velocity_data
    else:
        # 如果没有速度数据，创建全零的速度数据
        for _ in range(len(data_array)):
            ecef_velocity.append([0.0, 0.0, 0.0])
    
    # 提取姿态数据
    attitude_data = []
    has_attitude_data = attitude_cols and all(col >= 0 for col in attitude_cols)
    if has_attitude_data:
        for row in data_array:
            att_row = []
            for col in attitude_cols:
                if col < len(row):
                    att_row.append(row[col])
                else:
                    att_row.append(0.0)
            attitude_data.append(att_row)
    else:
        # 如果没有姿态数据，创建全零的姿态数据
        for _ in range(len(data_array)):
            attitude_data.append([0.0, 0.0, 0.0])
    
    # 合并所有数据，确保格式一致性
    output_data = []
    for i in range(len(data_array)):
        row_data = []
        
        # 添加时间戳（必须存在）
        if has_time_data:
            row_data.extend(time_data[i])
        else:
            row_data.extend([0.0, 0.0])  # 默认2列时间戳（周, 周内秒）
        
        # 添加位置数据
        if has_position_data:
            row_data.extend(ecef_position[i])
        else:
            row_data.extend([0.0, 0.0, 0.0])
        
        # 添加速度数据
        if has_velocity_data:
            row_data.extend(ecef_velocity[i])
        else:
            row_data.extend([0.0, 0.0, 0.0])
        
        # 添加姿态数据
        if has_attitude_data:
            row_data.extend(attitude_data[i])
        else:
            row_data.extend([0.0, 0.0, 0.0])
        
        output_data.append(row_data)
    
    return np.array(output_data)

def convert_reference_file():
    """将任意格式的参考文件转换为标准参考文件"""
    
    # 使用全局配置
    input_file = INPUT_FILE
    skip_lines = SKIP_LINES
    time_system = TIME_SYSTEM
    input_time_format = INPUT_TIME_FORMAT
    time_cols = TIME_COLS
    position_cols = POSITION_COLS
    position_format = POSITION_FORMAT
    velocity_cols = VELOCITY_COLS
    velocity_format = VELOCITY_FORMAT
    velocity_coefficients = VELOCITY_COEFFICIENTS
    attitude_cols = ATTITUDE_COLS
    output_file = OUTPUT_FILE
    verbose = VERBOSE
    
    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在")
        return False
    
    # 设置默认输出文件名
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}.truth"
    
    try:
        # 根据文件扩展名选择读取方法
        file_ext = os.path.splitext(input_file)[1].lower()
        
        if verbose:
            print(f"读取文件: {input_file}")
            print(f"文件类型: {file_ext}")
            print(f"时间系统: {time_system}")
            print(f"输入时间格式: {input_time_format}")
            print(f"输出时间格式: {OUTPUT_TIME_FORMAT}")
        
        if file_ext == '.csv':
            data_array = read_csv_file(input_file, skip_lines)
        elif file_ext in ['.xlsx', '.xls']:
            data_array = read_excel_file(input_file, skip_lines)
        else:
            # 默认为文本文件
            data_array = read_text_file(input_file, skip_lines)
        
        if data_array is None:
            print("错误: 无法读取文件数据")
            return False
        
        # 检查数据数组中的数据类型
        if verbose:
            print(f"数据数组形状: {data_array.shape}")
            print(f"数据数组数据类型: {data_array.dtype}")
        
        if verbose:
            print(f"读取数据成功，数据维度: {data_array.shape}")
        
        # 提取并转换数据
        output_data = extract_and_convert_data(
            data_array, time_cols, input_time_format, time_system,
            position_cols, position_format,
            velocity_cols, velocity_format, velocity_coefficients, attitude_cols
        )
        
        if output_data is None:
            print("错误: 数据提取和转换失败")
            return False
        
        if verbose:
            print(f"转换后数据维度: {output_data.shape}")
            print(f"时间系统: {time_system}")
            print(f"输入时间格式: {input_time_format}")
            print(f"时间戳列索引: {time_cols}")
            print(f"位置列索引: {position_cols} (格式: {position_format})")
            print(f"速度列索引: {velocity_cols} (格式: {velocity_format}, 系数: {velocity_coefficients})")
            print(f"姿态列索引: {attitude_cols}")
        
        # 生成右对齐的格式化字符串
        num_columns = output_data.shape[1] if len(output_data.shape) > 1 else 1
        fmt_string = generate_right_aligned_format(num_columns, FIELD_WIDTH, DECIMAL_PLACES)
        
        if verbose:
            print(f"生成的格式化字符串: '{fmt_string}'")
        
        # 保存为标准参考文件
        np.savetxt(output_file, output_data, fmt=fmt_string, delimiter=' ')
        
        print(f"\n转换完成!")
        print(f"输入文件: {input_file}")
        print(f"输出文件: {output_file}")
        print(f"时间系统: {time_system}")
        print(f"输入时间格式: {input_time_format}")
        print(f"输出时间格式: {OUTPUT_TIME_FORMAT}")
        print(f"跳过的行数: {skip_lines}")
        print(f"数据维度: {output_data.shape}")
        
        # 显示前几行数据预览
        if verbose and len(output_data) > 0:
            print(f"\n数据预览（前3行，右对齐格式）:")
            for i in range(min(3, len(output_data))):
                formatted_line = fmt_string % tuple(output_data[i])
                print(f"第{i+1}行: {formatted_line}")
        
        return True
        
    except Exception as e:
        print(f"转换过程中发生错误: {e}")
        return False

def main():
    """主函数"""
    print("开始转换参考文件...")
    print("=" * 50)
    
    # 显示当前配置
    print(f"输入文件: {INPUT_FILE}")
    print(f"跳过行数: {SKIP_LINES}")
    print(f"时间系统: {TIME_SYSTEM}")
    print(f"输入时间格式: {INPUT_TIME_FORMAT}")
    print(f"输出时间格式: {OUTPUT_TIME_FORMAT}")
    print(f"时间戳列索引: {TIME_COLS}")
    print(f"位置列索引: {POSITION_COLS} (格式: {POSITION_FORMAT})")
    print(f"速度列索引: {VELOCITY_COLS} (格式: {VELOCITY_FORMAT}, 系数: {VELOCITY_COEFFICIENTS})")
    print(f"姿态列索引: {ATTITUDE_COLS}")
    print(f"输出格式: 右对齐，字段宽度{FIELD_WIDTH}，小数位数{DECIMAL_PLACES}")
    print("=" * 50)
    
    # 执行转换
    success = convert_reference_file()
    
    if success:
        print("\n转换成功完成!")
    else:
        print("\n转换失败!")

if __name__ == "__main__":
    main()