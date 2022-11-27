from numpy import asarray, sqrt

is_nurbs = True

#
# First Square
poly_order_u = 1
poly_order_v = 1

u_knots = asarray([
    0,
    0,
    1,
    1
])

v_knots = asarray([
    0,
    0,
    1,
    1
])

control_points = asarray([
    [
        [1.,  0., 0., 1.],
        [1.,  0., 1., 1.]
    ],
    [
        [1.,  1., 0., 1.],
        [1.,  1., 1., 1.]
    ],
])

#
# Second Square
poly_order_u2 = 1
poly_order_v2 = 1

u_knots2 = asarray([
    0,
    0,
    1,
    1
])

v_knots2 = asarray([
    1,
    1,
    2,
    2
])

control_points2 = asarray([
    [
        [0.,  0., 1., 1.],
        [1.,  0., 1., 1.]
    ],
    [
        [0.,  1., 1., 1.],
        [1.,  1., 1., 1.]
    ],
])