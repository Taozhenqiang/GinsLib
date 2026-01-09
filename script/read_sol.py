# read_sol.py
import numpy as np

def read_solution(navfile, skip_lines=28, row=10000, col=18, sample=5e-2):
    """
    读取文件并返回处理后的数据。
    
    Parameters:
    - navfile: 输入文件路径和文件名。
    - skip_lines: 跳过文件开头的行数（默认为28行）。
    - row: 预设数组大小（默认为10000行）。
    - col: 文件的列数（默认为18列）。
    
    Returns:
    - data: 处理后的数据（NumPy 数组）。
    """
    try:
        with open(navfile, 'r') as file:
            # 跳过前skip_lines行
            for _ in range(skip_lines):
                file.readline()
            
            # 创建一个空的NumPy数组来存储数据
            data = np.zeros((row, col))
            dttol = sample / 1e3
            len_data = 0

            # 逐行读取文件数据
            for line in file:
                if line and (line.startswith('%') or line.startswith('#')):  # 跳过注释行
                    continue

                line = line.strip()  # 去除两端空白字符
                if not line:
                    continue

                # 使用空格作为分隔符分割字符串
                sline = line.split()

                # 判断第二个元素是否为整数（相当于检查时间戳）
                try:
                    if abs(float(sline[1])-round(float(sline[1]))) > (sample/2.0+dttol):
                        continue
                except ValueError:
                    continue

                len_data += 1
                if len_data > row:
                    row *= 2  # 扩展数组的大小
                    data = np.resize(data, (row, col))  # 扩展数组

                # 检查sline长度是否足够，避免IndexError, GNSS(24), GNSS/INS(39)
                if len(sline) == 24:
                    mode = 'GNSS'
                elif len(sline) == 46:
                    mode = 'GNSS/INS'
                else:
                    continue  # 如果行数据不足，跳过该行

                # 提取特定的数据列
                # GNSS: gps week, sow, pos[x/y/z], ratio, vel[x/y/z]
                if mode == 'GNSS':
                    mline = [sline[0], sline[1], 
                             sline[2], sline[3], sline[4], sline[14],
                             sline[15], sline[16], sline[17], 
                             0, 0, 0,
                             0, 0, 0, 0, 0, 0]  
                else:  # GNSS/INS: gps week, sow, pos[x/y/z], ratio, vel[x/y/z], att[pitch/roll/heading], bg[x/y/z], ba[x/y/z]
                    mline = [sline[0], sline[1], 
                             sline[2], sline[3], sline[4], sline[14],
                             sline[15], sline[16], sline[17], 
                             sline[24], sline[25], sline[26],
                         sline[33], sline[34], sline[35], sline[36], sline[37], sline[38]]

                # 将提取的数据存入NumPy数组
                try:
                    data[len_data - 1, :] = np.array(mline, dtype=float)
                except ValueError:
                    continue

            # 如果数据少于初始预设行数，删除多余的行
            if len_data < row:
                data = data[:len_data, :]

            return data

    except FileNotFoundError:
        raise FileNotFoundError(f"无法打开文件: {navfile}")
    except Exception as e:
        print(f"发生错误: {e}")
        return None