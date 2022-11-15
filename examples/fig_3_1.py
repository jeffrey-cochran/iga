from numpy import asarray, linspace

poly_order=3
knots = asarray([0,0,0,0,1,1,1,1])
control_points = asarray([
    [0.,0.],
    [0.25,1.],
    [0.75,3.],
    [1.,0.]
])