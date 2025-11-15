import numpy as np

def read_ref(filename, skip_lines=22, row=1000000, col=10, ins_flag=1, ins_interval=1/200, GNSS_interval=1, idx=None):
    """
    读取导航参考文件并将数据转换为NumPy数组。

    参数:
    - filename: str, 输入的文件路径。
    - skip_lines: int, 跳过文件开头的行数，默认22。
    - row: int, 初始分配的行数，默认1000000。
    - col: int, 数据的列数，默认10。
    - ins_flag: int, 指示是否为INS数据，1表示是INS，0表示其他类型数据，默认1。
    - interval: float, INS采样间隔，默认1/200。
    - idx: list, 需要读取的列索引列表，默认读取所有列。

    返回:
    - data: numpy.ndarray, 读取并处理后的数据。
    """
    if idx is None:
        idx = list(range(col))  # 默认读取col列
    
    try:
        with open(filename, 'r') as file:
            # 跳过前面的行
            for _ in range(skip_lines):
                next(file)
            
            data = np.zeros((row, col))
            len_data = 0
            
            for line in file:
                line = line.strip()

                # 检查当前行的首字母是否为字符
                if line and line[0].isalpha():  # 如果首字母是字母，则跳过该行
                    continue

                if not line: # 跳过空行
                    continue
                
                sline = line.split()
                time = float(sline[1])

                # 过滤条件
                if ((ins_flag and abs(time - round(time)) > 0.5 * ins_interval) or (len_data > 0 and abs(time - data[len_data - 1, 0]) < 0.5 * GNSS_interval)
                    or (not ins_flag and time % 1 != 0)):
                    continue
                
                len_data += 1
                if len_data > row:
                    row *= 2  # 扩展数组大小
                    data = np.resize(data, (row, col))
                
                mline = [float(sline[i]) for i in idx]  
                data[len_data - 1, :] = mline
            
            # 剔除多余的行
            if len_data < row:
                data = data[:len_data, :]
            
        return data

    except Exception as e:
        print(f"无法打开或读取文件: {e}")
        return None

