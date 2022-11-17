from numpy import asarray, sqrt

is_nurbs = True

poly_order_u = 1
poly_order_v = 2

u_knots = asarray([
    0,
    0,
    1,
    1
])

v_knots = asarray([
    0,
    0,
    0,
    .5,
    .5,
    1,
    1,
    1
])

a = sqrt(2.)/2.
control_points = asarray([
    [
        [1., 1., 0., 1.],
        [a, a, a, a],
        [1., 0., 1., 1.],
        [a, -a, a, a],
        [1., -1., 0, 1]
    ],
    [
        [1., 2., 0., 1.],
        [a, 2*a, 2*a, a],
        [1, 0., 2., 1.],
        [a, -2*a, 2*a, a],
        [1., -2., 0, 1]
    ]   
])