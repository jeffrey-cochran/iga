from numpy import asarray, sqrt

is_nurbs = True

poly_order_u = 1
poly_order_v = 1
poly_order_w = 2

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

w_knots = asarray([
    0,
    0,
    0,
    1,
    1,
    1
])

a = sqrt(2.)/2.
control_points = asarray([
    [
        [
            corner_1,
            [0.5, 0.5, 0, 1 ],
            [1, 0, 0, 1.]
        ],
        [
            [0,      2, 0, 1.],
            [2*a,  2*a, 0, a ],
            [2,      0, 0, 1.]
        ]
    ],
    [
        [
            [0, 1, 1, 1.],
            [0.5, 0.5, 1, 1 ],
            [1, 0, 1, 1.]
        ],
        [
            [0,     2, 1, 1.],
            [2*a, 2*a, a, a ],
            [2,     0, 1, 1.]
        ]
    ]
])