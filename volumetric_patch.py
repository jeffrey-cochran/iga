from patch import Patch
from surface_patch import SurfacePatch
from surface_union import SurfaceUnion
from numpy import typing, asarray, cumsum, cos, sin, dot, broadcast_to

class VolumetricPatch(Patch):

    def __init__(
        self, 
        planar_rotation=None,
        translation=None,
        object_name="obj_1",
        layer_name="lay_1",
        offset_vp=0,
        offset_v=0,
        offset_trim=0,
        *,
        knots_u,
        knots_v,
        knots_w,
        control_points,
        poly_order_u,
        poly_order_v,
        poly_order_w
    ):

        self.check_compatibility(
            knots_u,
            knots_v,
            knots_w,
            control_points,
            poly_order_u,
            poly_order_v,
            poly_order_w
        )

        poly_orders = [poly_order_u, poly_order_v, poly_order_w]
        knots = [knots_u, knots_v, knots_w]

        self.object_name = object_name
        self.layer_name = layer_name
        self.offset_vp = offset_vp
        self.offset_v = offset_v
        self.offset_trim = offset_trim

        self.keys = ['u','v','w']

        self.poly_order = {
            self.keys[i]: poly_orders[i] for i in range(3)
        }

        self.knot_vector = {
            self.keys[i]: knots[i] for i in range(3)
        }

        self.control_points = control_points

        if planar_rotation is not None:
            Rz = asarray([
                [cos(planar_rotation), -sin(planar_rotation), 0],
                [sin(planar_rotation),  cos(planar_rotation), 0],
                [                   0,                     0, 1]
            ])
            self.control_points[...,:3] = dot(Rz, self.control_points[:,:,:,:3].transpose(0,1,3,2)).transpose(1,2,3,0)

        if translation is not None:
            print(translation)
            weighted_translation = (
                broadcast_to(translation, self.control_points[...,:3].shape).transpose(3,0,1,2) * self.control_points[...,-1]
            ).transpose(1,2,3,0)
            self.control_points[...,:3] += weighted_translation

        self.bounding_surface = self.init_surfaces()
        (
            self.num_v,
            self.num_vp,
            self.num_trimming_curves
        ) = self.bounding_surface.object_counts()

        return

    def check_compatibility(
        self,
        knots_u,
        knots_v,
        knots_w,
        control_points,
        poly_order_u,
        poly_order_v,
        poly_order_w
    ):

        if len(control_points.shape) != 4:
            raise ValueError(f"The control point mesh must be 3D, but it's {len(control_points.shape)}D.")

        knot_dims = [len(z) for z in [knots_u, knots_v, knots_w]]
        poly_orders = [poly_order_u, poly_order_v, poly_order_w]
        for i in range(3):
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

    def set_offsets(self,*,
        offset_v, 
        offset_vp, 
        offset_trim
    ):
        self.offset_v = offset_v
        self.offset_vp = offset_vp
        self.offset_trim = offset_trim

        self.bounding_surface = self.init_surfaces()
        (
            self.num_v,
            self.num_vp,
            self.num_trimming_curves
        ) = self.bounding_surface.object_counts()
        return

    def dumps(self):
        return self.bounding_surface.dumps()

    def init_surfaces(self) -> SurfaceUnion:
        #
        # Direction of increasing u,v,w in the tensor product of the knot vectors
        #
        #            w+
        #           /
        #          +---v+
        #          |
        #          u+
        #
        #      +--------+
        #     /        /|
        #    /        / |
        #   +--------+  |
        #   |        |  |
        #   |        |  +
        #   |        | /
        #   |        |/
        #   +--------+
        #
        # Since u,v,w lie in [0,1], the surfaces correspond to those 
        # combinations for which one of u,v,w is fixed in {0,1}:
        #
        #   Uvw = 0vw = 1
        #   uVw = u0w = 2
        #   uvW = uv0 = 3
        #   uvW = uv1 = 4
        #   uVw = u1w = 5
        #   Uvw = 1vw = 6
        # 
        # I've numbered them this way so that the indices of two 
        # surfaces add to 7 if and only if they aren't adjacent (like a D6)
        #

        # self.surface_keys = [
        #     ('u', 0),
        #     ('v', 0),
        #     ('w', 0),
        #     ('w', 1),
        #     ('v', 1),
        #     ('u', 1)
        # ]

        return SurfaceUnion(
            object_name=self.object_name, 
            layer_name=self.layer_name,
            offset_trim=self.offset_trim,
            offset_v=self.offset_v,
            offset_vp=self.offset_vp,
            surface_patches=[
                SurfacePatch(
                    knots_u=self.knot_vector['v'], 
                    knots_v=self.knot_vector['w'], 
                    control_points=self.control_points[0,:,:,:], 
                    poly_order_u=self.poly_order['v'], 
                    poly_order_v=self.poly_order['w']
                ),
                SurfacePatch(
                    knots_u=self.knot_vector['u'], 
                    knots_v=self.knot_vector['w'], 
                    control_points=self.control_points[:,0,:,:], 
                    poly_order_u=self.poly_order['u'], 
                    poly_order_v=self.poly_order['w']
                ),
                SurfacePatch(
                    knots_u=self.knot_vector['u'], 
                    knots_v=self.knot_vector['v'], 
                    control_points=self.control_points[:,:,0,:], 
                    poly_order_u=self.poly_order['u'], 
                    poly_order_v=self.poly_order['v']
                ),
                SurfacePatch(
                    knots_u=self.knot_vector['u'], 
                    knots_v=self.knot_vector['v'], 
                    control_points=self.control_points[:,:,-1,:], 
                    poly_order_u=self.poly_order['u'], 
                    poly_order_v=self.poly_order['v']
                ),            
                SurfacePatch(
                    knots_u=self.knot_vector['u'], 
                    knots_v=self.knot_vector['w'], 
                    control_points=self.control_points[:,-1,:,:], 
                    poly_order_u=self.poly_order['u'], 
                    poly_order_v=self.poly_order['w']
                ),
                SurfacePatch(
                    knots_u=self.knot_vector['v'], 
                    knots_v=self.knot_vector['w'], 
                    control_points=self.control_points[-1,:,:,:], 
                    poly_order_u=self.poly_order['v'], 
                    poly_order_v=self.poly_order['w']
                )
            ]
        ) 