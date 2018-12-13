import numpy as np
import scipy.optimize
import functools
from scipy.spatial import distance


def atom_coordinates(atoms):
    return np.array([x.get_coord() for x in atoms])


def axis_sum(arr):
    values = [0] * arr.shape[1]

    for i in range(0, arr.shape[1]):
        values[i] = np.sum(arr[:, i])

    return np.array(values)


def point_in_line(p1, p2, d):
    v = p2 - p1
    u = v / np.linalg.norm(v)
    return p1 + d * u


def orthog_point(point, normal, t):
    return point + (normal * t)


def centroid(arr):
    return axis_sum(arr) / arr.shape[0]


def euclidean_distance(p1, p2):
    return distance.euclidean(p1, p2)


def angle(p1, p2):
    cosAngle = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    arcosAngle = np.arccos(np.clip(cosAngle, -1, 1))
    return np.degrees(arcosAngle)


def norm_vector(p1, p2):
    p1 = p1 / np.linalg.norm(p1)
    p2 = p2 / np.linalg.norm(p2)
    return p2 - p1


def to_quad1(angle):
    if (angle > 90 and angle <= 180):
        return 180 - angle
    elif (angle > 180 and angle <= 270):
        return angle - 180
    elif (angle > 270 and angle <= 360):
        return 360 - angle
    else:
        return angle


def plane(x, y, params):
    a = params[0]
    b = params[1]
    c = params[2]
    z = a * x + b * y + c
    return z


def error(params, points):
    result = 0
    for (x, y, z) in points:
        plane_z = plane(x, y, params)
        diff = abs(plane_z - z)
        result += diff**2
    return result


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]


def calc_normal(points):
    fun = functools.partial(error, points=points)
    params0 = [0, 0, 0]

    # It was found a Warning when executing the next line.
    # The Warning is:
    #       optimize.py:994: RuntimeWarning:
    #              divide by zero encountered in double_scalars
    #       rhok = 1.0 / (numpy.dot(yk, sk))
    #
    # However, according to my analysis, Numpy is able to deal
    # with such Warning and the result is generated properly.
    res = scipy.optimize.minimize(fun, params0)

    a = res.x[0]
    b = res.x[1]

    xs, ys, zs = zip(*points)

    normal = np.array(cross([1, 0, a], [0, 1, b]))

    return normal
