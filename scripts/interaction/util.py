import math
import numpy as np
from scipy.spatial import distance


def get_coordinatesnp(atoms):
    return np.array([x.get_coord() for x in atoms])


def get_sum_of_axis(arr):
    values = [0] * arr.shape[1]

    for i in range(0, arr.shape[1]):
        values[i] = np.sum(arr[:, i])

    return np.array(values)


def calc_centroidnp(arr):
    return get_sum_of_axis(arr) / arr.shape[0]


def calc_euclidean_distance(p1, p2):
    return distance.euclidean(p1, p2)


def norm_vector(p1, p2):
    p1 = p1 / np.linalg.norm(p1)
    p2 = p2 / np.linalg.norm(p2)
    return p2 - p1


def calc_angle(p1, p2):
    cosAngle = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    arcosAngle = np.arccos(np.clip(cosAngle, -1, 1))
    return np.degrees(arcosAngle)
    # return math.degrees(math.acos(cosAngle))


def to_quad1(angle):
    if (angle > 90 and angle <= 180):
        return 180 - angle
    elif (angle > 180 and angle <= 270):
        return angle - 180
    elif (angle > 270 and angle <= 360):
        return 360 - angle
    else:
        return angle


def get_point_in_line(p1, p2, d):
    v = p2 - p1
    u = v / np.linalg.norm(v)
    return p1 + d * u
