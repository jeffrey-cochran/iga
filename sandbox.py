import examples.ex_annulus as ex
from utils import surface_point, plot_bspline
from numpy import linspace
import matplotlib.pyplot as plt

parameter_values = linspace(0,1,100)
xyz = surface_point(
    parameter_values,
    ex.u_knots,
    ex.v_knots,
    ex.control_points,
    ex.poly_order_u,
    ex.poly_order_v,
    is_nurbs=True
)

X = xyz[:,:,0]
Y = xyz[:,:,1]
Z = xyz[:,:,2]

# print(z)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect([1,1,1])
surf = ax.plot_surface(X, Y, Z)

# plot_bspline(ex.knots, ex.control_points, ex.poly_order, is_nurbs=True )

plt.show()