import numpy as np

def cal_gravity(lat,h):
    """
    计算当地重力加速度 (m/s²) \n
    :param lat: 纬度 (rad)
    :param h: 海拔高度 (m)
    :return: 重力加速度 (m/s²) \n
    ref: I2NAV of WuHan University
    """

    a = 6378137.0           # 地球椭球长半轴 (m)
    b = 6356752.3141        # 地球椭球短半轴 (m)
    wie = 7.292115e-5       # 地球自转角速度 (rad/s)
    f = 1.0/298.257222101   # 地球扁率
    GM = 3.986005e14        # 地球引力常数 (m³/s²)
    gamma_a = 9.7803267715  # 赤道重力加速度 (m/s²)
    gamma_b = 9.8321863685  # 极点重力加速度 (m/s²)

    gamma = (a * gamma_a * np.cos(lat)**2 + b * gamma_b * np.sin(lat)**2) / np.sqrt(a**2*np.cos(lat)**2 + b**2*np.sin(lat)**2)
    m = wie**2 *a**2 * b / GM;

    g = gamma * (1 - 2 / a * (1 + f + m - 2 * f * np.sin(lat)**2) * h + 3 / a**2 * h**2)

    return g
