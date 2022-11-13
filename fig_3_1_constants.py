from numpy import asarray, linspace

knots = asarray([0,0,0,0,1,1,1,1])
parameter_values = linspace(0.0001,0.9999,1000)
control_points = asarray([
    [0.,0.],
    [0.25,1.],
    [0.75,3.],
    [1.,0.]
])