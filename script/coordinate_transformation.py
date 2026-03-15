import numpy as np

def ecef2pos(r):
    """
    Convert ECEF (Earth-Centered Earth-Fixed) coordinates to geographic coordinates (latitude, longitude, height).
    :param r: A numpy array or list with the ECEF coordinates [x, y, z].
    :return: pos: A numpy array containing [latitude, longitude, height].
    """
    FE_WGS84 = 1.0 / 298.257223563  # Flattening factor of WGS84
    RE_WGS84 = 6378137  # Earth's radius in meters (WGS84)

    e2 = FE_WGS84 * (2.0 - FE_WGS84)  # Square of eccentricity
    r2 = np.dot(r[0:2], r[0:2])  # r2 = x^2 + y^2
    v = RE_WGS84  # Initial approximation for the radius
    z = r[2]  # z-coordinate
    zk = 0  # Previous z value

    # Iterative solution for height and latitude
    while abs(z - zk) > 1E-4:
        zk = z
        sinp = z / np.sqrt(r2 + z * z)  # Sine of the latitude
        v = RE_WGS84 / np.sqrt(1.0 - e2 * sinp * sinp)  # Ellipsoid radius at latitude
        z = r[2] + v * e2 * sinp  # Adjusted z

    # Compute latitude, longitude, and height
    if r2 > 1E-12:
        pos = np.zeros(3)
        pos[0] = np.arctan(z / np.sqrt(r2))  # Latitude
        pos[1] = np.arctan2(r[1], r[0])  # Longitude
    else:
        if r[2] > 0:
            pos = np.array([np.pi / 2, 0, 0])  # North Pole
        else:
            pos = np.array([-np.pi / 2, 0, 0])  # South Pole

    pos[2] = np.sqrt(r2 + z * z) - v  # Height above the ellipsoid

    return pos

def xyz2enu(pos):
    """
    计算ECEF到局部坐标系的变换矩阵
    :param pos: 大地坐标数组 [纬度, 经度] (弧度)
    :return: ECEF到局部坐标变换矩阵 (3x3 numpy数组)
    """
    sinp = np.sin(pos[0])  # 纬度正弦
    cosp = np.cos(pos[0])  # 纬度余弦
    sinl = np.sin(pos[1])  # 经度正弦
    cosl = np.cos(pos[1])  # 经度余弦
    
    # 创建3x3变换矩阵（按列优先顺序存储）
    E = np.zeros((3, 3))
    
    # 第一列
    E[0, 0] = -sinl
    E[1, 0] = -sinp * cosl
    E[2, 0] = cosp * cosl
    
    # 第二列
    E[0, 1] = cosl
    E[1, 1] = -sinp * sinl
    E[2, 1] = cosp * sinl
    
    # 第三列
    E[0, 2] = 0.0
    E[1, 2] = cosp
    E[2, 2] = sinp
    
    return E

def xyz2blh(xyz):
    """
    Convert cartesian coordinates (XYZ) to geodetic coordinates (latitude, longitude, height).
    :param xyz: A numpy array or list with the ECEF coordinates [x, y, z].
    :return: A tuple (blh, Cne) where blh is a numpy array containing [latitude, longitude, height]
             and Cne is the transformation matrix (if requested).
    """
    # Convert ECEF (XYZ) to BLH using ecef2pos
    pos = ecef2pos(xyz)  # Assuming ecef2pos is defined elsewhere
    blh = pos  # Latitude, longitude, height

    B = pos[0]  # Latitude
    L = pos[1]  # Longitude
    sinB = np.sin(B)
    cosB = np.cos(B)
    sinL = np.sin(L)
    cosL = np.cos(L)

    Cne = None
    if len(blh) == 3:  # If transformation matrix is needed
        # Transformation matrix from ECEF (XYZ) frame to ENU frame
        Cne = np.array([
            [-sinL, cosL, 0],
            [-sinB * cosL, -sinB * sinL, cosB],
            [cosB * cosL, cosB * sinL, sinB]
        ])

    return blh, Cne

def pos2ecef(pos):
    """
    将大地坐标转换为ECEF坐标
    :param pos: 大地坐标数组 [纬度, 经度, 高度] (弧度, 米)
    :return: ECEF坐标数组 [x, y, z] (米)
    """
    # WGS84椭球参数
    FE_WGS84 = 1.0 / 298.257223563  # WGS84扁率因子
    RE_WGS84 = 6378137.0  # WGS84地球半径 (米)
    
    sinp = np.sin(pos[0]*np.pi/180.0)  # 纬度正弦
    cosp = np.cos(pos[0]*np.pi/180.0)  # 纬度余弦
    sinl = np.sin(pos[1]*np.pi/180.0)  # 经度正弦
    cosl = np.cos(pos[1]*np.pi/180.0)  # 经度余弦
    
    e2 = FE_WGS84 * (2.0 - FE_WGS84)  # 偏心率平方
    v = RE_WGS84 / np.sqrt(1.0 - e2 * sinp * sinp)  # 地球半径
    
    r = np.zeros(3)
    r[0] = (v + pos[2]) * cosp * cosl  # x坐标
    r[1] = (v + pos[2]) * cosp * sinl  # y坐标
    r[2] = (v * (1.0 - e2) + pos[2]) * sinp  # z坐标
    
    return r
