import re
import os
import argparse
from collections import Counter

# ========================= 配置区域 =========================
# 在这里修改配置参数（如果使用命令行参数，这些配置将被覆盖）

# 输入文件路径
INPUT_FILE = "./GNSS/LG69T_Vehicle_complex_20250414/LG69TAP01-TEMP0711-1HZ_INS.dat"

# 输出文件路径
OUTPUT_FILE = "./GNSS/LG69T_Vehicle_complex_20250414/imu.txt"  # None表示自动生成

# 目标指令列表（支持多个指令）
TARGET_COMMANDS = ["$PQTMRAWIMU"]  # 示例：提取$PQTMDRPVA指令
# TARGET_COMMANDS = ["$PQTMDRPVA", "$PQTMINSPVA", "$PQTMODOMSG"]  # 提取多个指令

# 文件编码（如果遇到编码错误可以修改）
FILE_ENCODING = "utf-8"  # 可选: "gbk", "gb2312", "utf-8"

# 是否显示详细信息
VERBOSE = True

# 是否在输出中包含原始行（包括指令头）
INCLUDE_COMMAND_HEADER = True

# ========================= 函数定义 =========================

def parse_command_line(line, command_pattern):
    """
    使用正则表达式解析指令行
    """
    # 构建正则表达式模式
    pattern = rf"^{re.escape(command_pattern)}(.*)$"
    match = re.match(pattern, line.strip())
    
    if match:
        # 提取指令后面的数据部分
        data_part = match.group(1)
        
        # 使用多种分隔符分割数据
        # 支持逗号、星号、分号、空格等分隔符
        separators = r'[,\*;\s]+'
        elements = re.split(separators, data_part)
        
        # 过滤空字符串
        elements = [elem for elem in elements if elem]
        
        return elements
    
    return None

def extract_commands_from_file():
    """
    从文件中提取指定指令
    """
    # 使用全局配置
    input_file = INPUT_FILE
    output_file = OUTPUT_FILE
    target_commands = TARGET_COMMANDS
    encoding = FILE_ENCODING
    verbose = VERBOSE
    include_header = INCLUDE_COMMAND_HEADER
    
    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 不存在")
        return False
    
    # 设置默认输出文件名
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        command_names = "_".join([cmd.replace("$", "") for cmd in target_commands])
        output_file = f"{base_name}_{command_names}_extracted.txt"
    
    try:
        # 读取文件
        extracted_data = {}
        total_lines = 0

        # 读取文件内容，使用replace编码防止丢行
        with open(input_file, 'r', encoding=encoding, errors='replace') as f:
            for line_num, line in enumerate(f, 1):
                total_lines += 1
                line = line.strip()
                
                line = re.sub(r'^[^$]*', '', line)  # 删掉开头所有非 $ 字符

                # 跳过空行
                if not line:
                    continue
                
                # 检查是否包含目标指令
                for command in target_commands:
                    # 注意：这里使用包含关系检查，而不是以指令开头检查
                    if command in line:
                        # 解析指令行
                        elements = parse_command_line(line, command)
                        
                        if elements:
                            if command not in extracted_data:
                                extracted_data[command] = []
                            
                            if include_header:
                                # 包含指令头的完整数据
                                extracted_data[command].append([command] + elements)
                            else:
                                # 只包含数据部分
                                extracted_data[command].append(elements)
        
        if verbose:
            print(f"\n处理完成!")
            print(f"总行数: {total_lines}")
        
        # 保存提取结果
        with open(output_file, 'w', encoding='utf-8') as f_out:
            for command, data_list in extracted_data.items():
                for data in data_list:
                    # 将数据转换为字符串格式
                    if include_header:
                        # 重新构建原始格式（指令头 + 数据）
                        line_content = ",".join(data)
                    else:
                        # 只保存数据部分
                        line_content = ",".join(data)
                    
                    f_out.write(line_content + "\n")
        
        # 显示提取统计
        for command, data in extracted_data.items():
            print(f"{command}: 提取 {len(data)} 行")
            
            # 显示前3行示例
            if verbose and data:
                print(f"  前3行示例:")
                for i in range(min(3, len(data))):
                    print(f"    第{i+1}行: {data[i][:5]}...")  # 只显示前5个元素
        
        return True
        
    except Exception as e:
        print(f"提取过程中发生错误: {e}")
        return False

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='从文本文件中提取指定指令')
    parser.add_argument('-i', '--input', type=str, help='输入文件路径')
    parser.add_argument('-o', '--output', type=str, help='输出文件路径')
    parser.add_argument('-c', '--commands', type=str, help='目标指令(逗号分隔，如:$PQTMDRPVA,$PQTMINSPVA)')
    parser.add_argument('-e', '--encoding', type=str, help='文件编码')
    parser.add_argument('--no-header', action='store_true', help='不包含指令头')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细信息')
    
    return parser.parse_args()

def main():
    """主函数"""
    # 解析命令行参数
    args = parse_arguments()
    
    # 更新配置参数
    global INPUT_FILE, OUTPUT_FILE, TARGET_COMMANDS, FILE_ENCODING, VERBOSE, INCLUDE_COMMAND_HEADER
    
    if args.input:
        INPUT_FILE = args.input
    if args.output:
        OUTPUT_FILE = args.output
    if args.commands:
        TARGET_COMMANDS = [cmd.strip() for cmd in args.commands.split(',')]
    if args.encoding:
        FILE_ENCODING = args.encoding
    if args.no_header:
        INCLUDE_COMMAND_HEADER = False
    if args.verbose:
        VERBOSE = args.verbose
    
    print("开始提取指定指令...")
    print("=" * 50)
    
    # 显示当前配置
    print(f"输入文件: {INPUT_FILE}")
    print(f"输出文件: {OUTPUT_FILE}")
    print(f"目标指令: {TARGET_COMMANDS}")
    print(f"文件编码: {FILE_ENCODING}")
    print(f"包含指令头: {INCLUDE_COMMAND_HEADER}")
    print("=" * 50)
    
    # 执行提取
    success = extract_commands_from_file()
    
    if success:
        print("\n指令提取成功!")
    else:
        print("\n指令提取失败!")

if __name__ == "__main__":
    main()