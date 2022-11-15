from numpy import asarray, linspace, sqrt

is_nurbs = True
poly_order = 2

knots = asarray(
    [
        0,
        0,
        0,
        1,
        1,
        1
])

a = sqrt(2.)/2.
weights = asarray([
    1.,
    a,
    1.
])

uw_control_points = asarray([
    [0., 1., 1.],
    [1., 1., 1. ],
    [1., 0., 1.]
])
control_points = (uw_control_points.transpose()*weights).transpose()