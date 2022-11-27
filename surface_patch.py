from patch import Patch
from trim_curve import TrimCurve
from numpy import asarray
import tempfile

class SurfacePatch(Patch):

    def __init__(
        self, 
        offset_vp=0,
        offset_v=0,
        offset_trim=0,
        *,
        knots_u,
        knots_v,
        control_points,
        poly_order_u,
        poly_order_v

    ):

        self.check_compatibility(
            knots_u,
            knots_v,
            control_points,
            poly_order_u,
            poly_order_v
        )

        poly_orders = [poly_order_u, poly_order_v]
        knots = [knots_u, knots_v]

        self.offset_vp = offset_vp
        self.offset_v = offset_v
        self.offset_trim = offset_trim

        self.keys = ['u','v']

        self.poly_order = {
            self.keys[i]: poly_orders[i] for i in range(2)
        }

        self.knot_vector = {
            self.keys[i]: knots[i] for i in range(2)
        }

        self.control_points = control_points

        self.num_trim_curves = 4
        self.num_v = self.control_points[...,0].size
        self.num_vp = self.num_trim_curves*2

        self.trimming_curves = self.init_trim_curves()

        return

    def check_compatibility(
        self,
        knots_u,
        knots_v,
        control_points,
        poly_order_u,
        poly_order_v
    ):

        if len(control_points.shape) != 3:
            raise ValueError(f"The control point mesh must be 3D, but it's {len(control_points.shape)}D.")

        knot_dims = [len(z) for z in [knots_u, knots_v]]
        poly_orders = [poly_order_u, poly_order_v]
        for i in range(2):
            self.check_compatibility_1d(
                knot_dims[i],
                control_points.shape[i],
                poly_orders[i],
                i
            )

        return

    def check_compatibility_1d(
        self,
        knots_dim,
        control_points_dim,
        poly_order,
        direction
    ):
        correct_dimension = knots_dim - poly_order - 1
        if control_points_dim != correct_dimension:
            raise ValueError(
                f"The control mesh should have {correct_dimension} control points in the {direction}-direction, but it has {control_points_dim} instead."
            )
            
        return

    def dumps(self, include_trim_loop=True):

        trim_loop_curves_string = "" 
        trim_loop_statement_string = ""

        if include_trim_loop:
            trim_loop_statement_string = "trim "
            temp = []
            trim_offset = 1
            for trimming_curve in self.trimming_curves:
                temp += [
                    str(trimming_curve.u0), 
                    str(trimming_curve.u1), 
                    str(self.offset_trim + trim_offset)
                ]
                trim_loop_curves_string += trimming_curve.dumps()
                trim_offset += 1
            temp_str = " ".join(temp)
            trim_loop_statement_string += temp_str


        knots_u = [str(a) for a in self.knot_vector['u'].tolist()]
        knots_v = [str(b) for b in self.knot_vector['v'].tolist()]
        knot_bounds = " ".join([
            str(self.knot_vector['u'][0]),
            str(self.knot_vector['u'][1]),
            str(self.knot_vector['v'][0]),
            str(self.knot_vector['v'][1])
        ])
        knot_string_u = " ".join(knots_u)
        knot_string_v = " ".join(knots_v)
        control_point_string = ""
        for j in range(self.control_points.shape[1]):
            for i in range(self.control_points.shape[0]):
                cp = self.control_points[i,j,:].copy()
                print("#-------------------------------------")
                print(cp)
                cp[:3] = cp[:3] / cp[3]
                print(cp)
                control_point_string += f"v {cp[0]} {cp[1]} {cp[2]} {cp[3]}\n"
        idx = [
            str(i+self.offset_v+1) for i in range(
                    self.control_points.shape[0] * self.control_points.shape[1]
            )
        ]
        vertex_string = " ".join(idx)

        return f"""# Begin surface dump

{"# Begin trim loop" if include_trim_loop else ""}
{trim_loop_curves_string}

{control_point_string}
cstype rat bspline
deg {self.poly_order['u']} {self.poly_order['v']}
surf {knot_bounds} {vertex_string}
parm u {knot_string_u}
parm v {knot_string_v}
{trim_loop_statement_string}
end
"""

    def set_offsets(self,*,
        offset_v, 
        offset_vp, 
        offset_trim
    ):
        self.offset_v = offset_v
        self.offset_vp = offset_vp
        self.offset_trim = offset_trim

        self.trimming_curves = self.init_trim_curves()
        return

    def init_trim_curves(self):

        u_min = self.knot_vector['u'].min()
        u_max = self.knot_vector['u'].max()

        v_min = self.knot_vector['v'].min()
        v_max = self.knot_vector['v'].max()

        trim_loop = [
            TrimCurve(
                point_1=(u_min, v_min),
                point_2=(u_max, v_min),
                degree=1,
                knots=asarray([u_min,u_min,u_max,u_max]),
                offset=self.offset_vp
            ),
            TrimCurve(
                point_1=(u_max, v_min),
                point_2=(u_max, v_max),
                degree=1,
                knots=asarray([v_min,v_min,v_max,v_max]),
                offset=self.offset_vp+2
            ),
            TrimCurve(
                point_1=(u_max, v_max),
                point_2=(u_min, v_max),
                degree=1,
                knots=asarray([u_min,u_min,u_max,u_max]),
                offset=self.offset_vp+4
            ),
            TrimCurve(
                point_1=(u_min, v_max),
                point_2=(u_min, v_min),
                degree=1,
                knots=asarray([v_min,v_min,v_max,v_max]),
                offset=self.offset_vp+6
            )
        ]
        return trim_loop
