from fig_3_1_constants import knots, parameter_values, control_points
from utils import curve_point
import matplotlib.pyplot as pl

# print(parameter_values.shape)
z = curve_point(parameter_values, knots, control_points, 3)
x = z[:,0]
y = z[:,1]

pl.plot(x,y)
pl.show()
