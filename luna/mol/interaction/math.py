import numpy as np
import scipy.optimize
import functools
from scipy.spatial import distance


def atom_coordinates(atoms):
    return np.array([x.get_coord() for x in atoms])


def ob_atom_coordinates(ob_atoms):
    return np.array([(a.GetX(), a.GetY(), a.GetZ()) for a in ob_atoms])


def axis_sum(arr, decimals=3):
    values = [0] * arr.shape[1]
    for i in range(0, arr.shape[1]):
        values[i] = np.around(np.sum(arr[:, i]), decimals)
    return np.array(values)


def point_in_line(p1, p2, d):
    v = p2 - p1
    u = v / np.linalg.norm(v)
    return p1 + d * u


def orthog_point(point, normal, t):
    return point + (normal * t)


def centroid(arr, decimals=3):
    return np.around(axis_sum(arr) / arr.shape[0], decimals)


def euclidean_distance(p1, p2, decimals=3):
    return round(distance.euclidean(p1, p2), decimals)


def angle(p1, p2, decimals=3):
    normal_prod = np.linalg.norm(p1) * np.linalg.norm(p2)

    if normal_prod == 0:
        return np.float32(0)

    dot_prod = np.dot(p1, p2)
    cos_angle = dot_prod / normal_prod
    arcos_angle = np.arccos(np.clip(cos_angle, -1, 1))

    return np.around(np.degrees(arcos_angle), 3)


def norm_vector(p1, p2):
    p1 = p1 / np.linalg.norm(p1)
    p2 = p2 / np.linalg.norm(p2)
    return p2 - p1


def to_quad1(angle, decimals=3):
    if (angle > 90 and angle <= 180):
        return np.around((180 - angle), decimals)
    elif (angle > 180 and angle <= 270):
        return np.around((angle - 180), decimals)
    elif (angle > 270 and angle <= 360):
        return np.around((360 - angle), decimals)
    else:
        return np.around(angle, decimals)


def plane(x, y, params):
    a = params[0]
    b = params[1]
    c = params[2]

    return a * x + b * y + c


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


def calc_normal(points, decimals=3):
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

    return np.around(np.array(cross([1, 0, a], [0, 1, b])), decimals)
