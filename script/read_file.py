import numpy as np
from coordinate_transformation import pos2ecef, ecef2pos, xyz2enu

def read_pos(navfile, skip_lines=28, row=10000, col=18, sample=5e-2, sol_type='GINSLIB', sol_idx=[], pos_type='ECEF', vel_type='ECEF'):
    """
    读取文件并返回处理后的数据。
    
    Parameters:
    - navfile: 输入文件路径和文件名
    - skip_lines: 跳过文件开头的行数
    - row: 预设数组大小（默认为10000行）
    - col: 文件的列数（默认为18列）
    
    Returns:
    - data: 处理后的数据（NumPy 数组）
    """
    try:
        with open(navfile, 'r') as file:
            # 跳过前skip_lines行
            for _ in range(skip_lines):
                file.readline()
            
            # 创建一个空的NumPy数组来存储数据
            if sol_idx and sol_type != 'GINSLIB':
                col = len(sol_idx)
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
                if sol_type == 'GINSLIB':
                    if len(sline) > 5 and len(sline) < 47:
                        mode = 'GNSS'
                    elif len(sline) == 47 or len(sline) == 46:
                        mode = 'GNSS/INS'
                    else:
                        continue  # 如果行数据不足，跳过该行

                    # 提取特定的数据列
                    if mode == 'GNSS':
                        # GNSS: gps week, sow, pos[x/y/z]
                        if len(sline) < 18:
                            mline = [sline[0], sline[1], 
                                    sline[2], sline[3], sline[4], 
                                    0, 0, 0, 
                                    0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0] 
                        else:
                            # GNSS: gps week, sow, pos[x/y/z], vel[x/y/z], 
                            mline = [sline[0], sline[1], 
                                    sline[2], sline[3], sline[4], 
                                    sline[15], sline[16], sline[17], 
                                    0, 0, 0, 
                                    sline[5], sline[6], sline[14], sline[24],
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0]  
                    else:  # GNSS/INS: gps week, sow, pos[x/y/z], vel[x/y/z], att[pitch/roll/heading], 
                           # Q, ns, ratio, PDOP
                           # bg[x/y/z], ba[x/y/z]
                           # pos_std, vel_std, att_std, bg_std, ba_std
                        mline = [sline[0], sline[1], 
                                sline[2], sline[3], sline[4], 
                                sline[15], sline[16], sline[17], 
                                sline[25], sline[26], sline[27], 
                                sline[5], sline[6], sline[14], sline[24],
                                sline[34], sline[35], sline[36], sline[37], sline[38], sline[39],
                                sline[7], sline[8], sline[9], sline[18], sline[19], sline[20], sline[28], sline[29], sline[30],
                                sline[40], sline[41], sline[42], sline[43], sline[44], sline[45]]
                
                else:
                    # 正确提取指定索引的数据
                    mline = [sline[i] for i in sol_idx]

                    if pos_type == 'LLH':
                        # 转换为弧度并计算ECEF坐标
                        lat_rad = np.deg2rad(float(mline[2]))
                        lon_rad = np.deg2rad(float(mline[3]))
                        height = float(mline[4])
                        llh_pos = [lat_rad, lon_rad, height]
                        ecef_pos = pos2ecef(llh_pos)
                        # 更新位置数据
                        mline[2:5] = ecef_pos
                    elif pos_type == 'ECEF':
                        llh_pos = ecef2pos([float(x) for x in mline[2:5]])

                    if vel_type == 'ENU' or vel_type == 'NED':
                        # 获取位置坐标用于计算变换矩阵
                        Cne = xyz2enu(llh_pos)
                        
                        if vel_type == 'NED':
                            # NED到ENU转换: N->E, E->N, D->-U
                            mline[5:8] = [float(mline[6]), float(mline[5]), -float(mline[7])]
                        
                        # ENU到ECEF转换
                        vel_enu = np.array([float(x) for x in mline[5:8]])
                        vel_ecef = np.dot(Cne.T, vel_enu)
                        mline[5:8] = vel_ecef.tolist()
                    
                    # 补全数据列(11列, gps week, sow, pos[x/y/z], vel[x/y/z], att[pitch/roll/heading])
                    if len(mline) < 11:
                        for i in range(len(mline), 11):
                            mline.append(0.0)
                                       
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