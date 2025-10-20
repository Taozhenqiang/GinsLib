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
