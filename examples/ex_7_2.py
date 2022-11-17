from numpy import asarray, linspace, sqrt

is_nurbs = True
poly_order = 2

knots = asarray(
    [
        0,
        0,
        0,
        1/4.,
        1/4.,
        1/2.,
        1/2.,
        3/4.,
        3/4.,
        1,
        1,
        1
])

a = sqrt(2.)/2.
control_points = asarray([
        [ 2,  0, 1],
        [ 2*a,  2*a, a],
        [ 0,  2, 1],
        [-2*a,  2*a, a],
        [-2,  0, 1],
        [-2*a, -2*a, a],
        [ 0, -2, 1],
        [ 2*a, -2*a, a],
        [ 2,  0, 1]
])