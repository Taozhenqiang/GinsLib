import numpy as np
from ecef2pos import ecef2pos

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
