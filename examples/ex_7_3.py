from numpy import asarray, linspace, cos, pi

is_nurbs = True
poly_order = 2

knots = asarray([0,0,0,1./3.,1./3.,2./3.,2./3.,1,1,1])

a = cos(pi/6.)
weights = asarray([
    1.0,
    0.5,
    1.0,
    0.5,
    1.0,
    0.5,
    1.0
])

uw_control_points = asarray([
    [   a,   0.5, 1.],
    [  0.,   2.0, 1.],
    [  -a,   0.5, 1.],
    [-2*a,   -1., 1.],
    [   0,   -1., 1.],
    [ 2*a,   -1., 1.],
    [   a,   0.5, 1.]
])
control_points = (uw_control_points.transpose()*weights).transpose()