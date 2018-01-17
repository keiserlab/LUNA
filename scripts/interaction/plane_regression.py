import numpy as np
import scipy.optimize
import functools


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


def get_orthog_point(point, normal, t):
    return point + (normal * t)
