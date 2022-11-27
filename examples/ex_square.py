from numpy import asarray, sqrt

is_nurbs = True

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