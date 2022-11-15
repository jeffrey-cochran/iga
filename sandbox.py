import examples.ex_7_3 as ex
from utils import plot_bspline

plot_bspline(ex.knots, ex.control_points, ex.poly_order, is_nurbs=ex.is_nurbs, show_control_points=True)