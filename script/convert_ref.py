import numpy as np
import os
import re
import pandas as pd
import csv

# ========================= 配置区域 =========================
# 在这里修改配置参数

# 输入文件路径
# INPUT_FILE = "./EG320N_Vehicle_open_20250418_campus2/TC_2GNSS.txt"  
# INPUT_FILE = "./GNSS/LG69T_GNSS_open_20260211/ref.xlsx"
INPUT_FILE = "./GNSS/LG69T_GNSS_open_20260211/ref.csv"  

# 要跳过的行数（通常是表头或注释行）
SKIP_LINES = 1

# 要提取的数据列索引（从0开始计数，不包含非数字项）
# week, sec, ecef_pos[x/y/z], ecef_vel[x/y/z], att[pitch/roll/heading]
# tzq
EXTRACT_COLS = [0, 1, 2, 3, 4]
# EXTRACT_COLS = [0, 1, 2, 3, 4, 8, 9, 10, 14, 15, 16]
# EXTRACT_COLS = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 8]
# zf
# EXTRACT_COLS = [0, 1, 6, 7, 8, 12, 13, 14, 19, 20, 18]

# 输出文件路径（如果为None，则自动生成）
OUTPUT_FILE = "./GNSS/LG69T_GNSS_open_20260211/truth.truth"  # 例如: "output.truth"

# 是否显示详细信息
VERBOSE = True

# 输出格式配置
# 右对齐格式：设置字段宽度，数字右对齐
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
    
    # 返回出现次数最多的分隔符
    if space_count > comma_count and space_count > tab_count:
        return 'space'
    elif comma_count > space_count and comma_count > tab_count:
        return 'comma'
    elif tab_count > space_count and tab_count > comma_count:
        return 'tab'
    else:
        return 'space'  # 默认使用空格

def parse_line(line, delimiter):
    """
    根据分隔符解析行数据
    """
    if delimiter == 'space':
        return re.split(r'\s+', line.strip())
    elif delimiter == 'comma':
        return line.strip().split(',')
    elif delimiter == 'tab':
        return line.strip().split('\t')
    else:
        return line.strip().split()

def generate_right_aligned_format(num_columns, field_width=15, decimal_places=6):
    """
    生成右对齐的格式化字符串
    
    参数:
    - num_columns: 列数
    - field_width: 字段宽度
    - decimal_places: 小数位数
    
    返回:
    - 格式化字符串，如: '%15.6f %15.6f %15.6f'
    """
    # 生成每个字段的格式：%[宽度].[精度]f
    format_strings = [f"%{field_width}.{decimal_places}f" for _ in range(num_columns)]
    # 用空格连接所有格式
    return ' '.join(format_strings)

def read_csv_file(file_path, skip_lines=0, extract_cols=None):
    """
    读取CSV文件
    """
    try:
        # 使用pandas读取CSV文件
        df = pd.read_csv(file_path, skiprows=skip_lines, header=None)
        
        # 提取指定列
        if extract_cols:
            df = df.iloc[:, extract_cols]
        
        # 转换为numpy数组
        data_array = df.to_numpy()
        
        return data_array
    except Exception as e:
        print(f"读取CSV文件错误: {e}")
        return None

def read_excel_file(file_path, skip_lines=0, extract_cols=None, sheet_name=0):
    """
    读取Excel文件
    """
    try:
        # 使用pandas读取Excel文件
        df = pd.read_excel(file_path, skiprows=skip_lines, header=None, sheet_name=sheet_name)
        
        # 提取指定列
        if extract_cols:
            df = df.iloc[:, extract_cols]
        
        # 转换为numpy数组
        data_array = df.to_numpy()
        
        return data_array
    except Exception as e:
        print(f"读取Excel文件错误: {e}")
        return None

def read_text_file(file_path, skip_lines=0, extract_cols=None):
    """
    读取文本文件（原有功能）
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as f_in:
            lines = f_in.readlines()
        
        # 跳过指定行数
        if skip_lines >= len(lines):
            print(f"警告: 跳过的行数({skip_lines})超过文件总行数({len(lines)})")
            skip_lines = 0
        
        data_lines = lines[skip_lines:]
        
        # 检测分隔符（使用第一个非空行）
        delimiter = 'space'
        for line in data_lines:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('%'):  # 跳过注释行
                delimiter = detect_delimiter(line)
                break
        
        # 处理数据
        processed_data = []
        
        for i, line in enumerate(data_lines, skip_lines + 1):
            line = line.strip()
            
            # 跳过空行和注释行
            if not line or line.startswith('#') or line.startswith('%'):
                continue
            
            # 解析行数据
            try:
                elements = parse_line(line, delimiter)
                
                # 过滤掉非数字元素
                numeric_elements = []
                for elem in elements:
                    try:
                        numeric_elements.append(float(elem))
                    except ValueError:
                        continue  # 跳过非数字元素
                
                # 检查是否有足够的数字列
                if len(numeric_elements) >= max(extract_cols) + 1:
                    # 提取指定列
                    extracted_data = [numeric_elements[col] for col in extract_cols if col < len(numeric_elements)]
                    processed_data.append(extracted_data)
                
            except Exception as e:
                continue
        
        if not processed_data:
            return None
        
        # 转换为numpy数组
        data_array = np.array(processed_data)
        
        return data_array
    except Exception as e:
        print(f"读取文本文件错误: {e}")
        return None

def convert_reference_file():
    """
    将任意格式的参考文件转换为标准参考文件
    """
    
    # 使用全局配置
    input_file = INPUT_FILE
    skip_lines = SKIP_LINES
    extract_cols = EXTRACT_COLS
    output_file = OUTPUT_FILE
    verbose = VERBOSE
    field_width = FIELD_WIDTH
    decimal_places = DECIMAL_PLACES
    
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
        
        if file_ext == '.csv':
            data_array = read_csv_file(input_file, skip_lines, extract_cols)
        elif file_ext in ['.xlsx', '.xls']:
            data_array = read_excel_file(input_file, skip_lines, extract_cols)
        else:
            # 默认为文本文件
            data_array = read_text_file(input_file, skip_lines, extract_cols)
        
        if data_array is None:
            print("错误: 无法读取文件数据")
            return False
        
        if verbose:
            print(f"读取数据成功，数据维度: {data_array.shape}")
            print(f"提取列索引: {extract_cols}")
        
        # 生成右对齐的格式化字符串
        num_columns = data_array.shape[1] if len(data_array.shape) > 1 else 1
        fmt_string = generate_right_aligned_format(num_columns, field_width, decimal_places)
        
        if verbose:
            print(f"生成的格式化字符串: '{fmt_string}'")
        
        # 保存为标准参考文件（右对齐格式）
        np.savetxt(output_file, data_array, fmt=fmt_string, delimiter=' ')
        
        print(f"\n转换完成!")
        print(f"输入文件: {input_file}")
        print(f"输出文件: {output_file}")
        print(f"跳过的行数: {skip_lines}")
        print(f"提取的列索引: {extract_cols}")
        print(f"数据维度: {data_array.shape}")
        print(f"输出格式: 右对齐，字段宽度{FIELD_WIDTH}，小数位数{DECIMAL_PLACES}")
        
        # 显示前几行数据预览
        if verbose and len(data_array) > 0:
            print(f"\n数据预览（前3行，右对齐格式）:")
            for i in range(min(3, len(data_array))):
                # 使用相同的格式化字符串显示预览
                formatted_line = fmt_string % tuple(data_array[i])
                print(f"第{i+1}行: {formatted_line}")
        
        return True
        
    except Exception as e:
        print(f"转换过程中发生错误: {e}")
        return False

def main():
    """
    主函数
    """
    print("开始转换参考文件...")
    print("=" * 50)
    
    # 显示当前配置
    print(f"输入文件: {INPUT_FILE}")
    print(f"跳过行数: {SKIP_LINES}")
    print(f"提取列索引: {EXTRACT_COLS}")
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