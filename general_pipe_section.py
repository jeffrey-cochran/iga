from volumetric_patch import VolumetricPatch
from numpy import cos, tan, asarray
from numpy.linalg import norm
from enums import FlatFace

class GeneralPipeSection(VolumetricPatch):

    def __init__(self, 
        planar_rotation=None,
        translation=None,
        *,
        corners,
        flat_face
    ):

        poly_order_u = 1
        poly_order_v = 1
        poly_order_w = 2

        knot_vector_u = asarray([0, 0, 1, 1])
        knot_vector_v = asarray([0, 0, 1, 1])
        knot_vector_w = asarray([0, 0, 0, 1, 1, 1])

        control_points = self.compute_control_points(corners, flat_face)
        print(translation)
        super().__init__(
            knots_u=knot_vector_u,
            knots_v=knot_vector_v,
            knots_w=knot_vector_w,
            control_points=control_points,
            poly_order_u=poly_order_u,
            poly_order_v=poly_order_v,
            poly_order_w=poly_order_w,
            planar_rotation=planar_rotation,
            translation=translation
        )

        return

    def compute_control_points(self,
        corners,
        flat_face
    ):
        """
        Looking in the positive r-direction...
        
                    z
                    |  r
                    | /
                    |/
         theta------+
        
        ...the corners are provided in the following
        order.
        
        First, the interior corners (inner surface):
        
           /       /
          /       /
         1---a---4
         | /     | /
         |/      |/
         2---b---3
        
        and then the exterior corners (outer surface):
        
          5---c---8
         /|      /|
        / |     / |
          6---d---7
         /       /
        /       /

        The control points a, b, c, and d are computed  
        based on the angle and flat_face.

        --------------------------

        If flat_face == FlatFlace.NONE, it means this is a 
        regular pipe section, and both the surface 1-2-3-4
        and the surface 5-6-7-8 are *not* flat, but cylindrical
        sections.

        If flat_face == FlatFlace.INT, it means the surface
        1-2-3-4 is flat, but the surface 5-6-7-8 is a cylindrical
        section.

        If flat_face == FlatFlace.EXT, it means the surface
        1-2-3-4 is a cylindrical section, but the surface 5-6-7-8 
        is flat.
        """

        corner_1 = corners[0]
        corner_2 = corners[1]
        corner_3 = corners[2]
        corner_4 = corners[3]
        corner_5 = corners[4]
        corner_6 = corners[5]
        corner_7 = corners[6]
        corner_8 = corners[7]

        #
        # Default case is a regular pipe section, i.e., 
        # the only qaudrilateral surfaces are in the r-z plane
        #
        # Compute internal and external radii
        r_int = norm(corner_4[:2])
        r_ext = norm(corner_8[:2])
        #
        # Compute the mid points
        mid_point_a = (corner_1[:2] + corner_4[:2]) / 2.
        mid_point_b = (corner_2[:2] + corner_3[:2]) / 2.
        mid_point_c = (corner_5[:2] + corner_8[:2]) / 2.
        mid_point_d = (corner_6[:2] + corner_7[:2]) / 2.
        #
        # Compute the radii of the control points
        # NOTE: cos(theta/2) = ||mid_point_a|| / r_int = r_int / ||a||, 
        # since both o-4-a and o-mid_a-4 are right angles. This means that
        # ||a|| = r_int * (r_int / ||mid__point_a||). Similar reasoning holds 
        # for b, c, and d. 
        #                                               
        #       a                                        c
        #       |                                        |
        # 1---mid_a---4                            5---mid_c---8
        #  \    |    /                              \    |    /
        #   \   |   /                                \   |   /
        #    \  |  / <~~~~ ||4|| = r_int              \  |  /  <~~~~ ||8|| = r_ext
        #     \ | /                                    \ | /
        #      \|/                                      \|/
        #       o                                        o
        #
        # NOTE: we require that ||a|| = ||b|| and that ||c|| = ||d||
        # NOTE: all of these calculations are planar calculations in the x-y plane
        norm_a = r_int * (r_int / norm(mid_point_a))
        norm_b = norm_a
        norm_c = r_ext * (r_ext / norm(mid_point_c))
        norm_d = norm_c
        #
        # Compute the unweighted points a, b, c, d
        point_a = norm_a * (mid_point_a / norm(mid_point_a))
        point_b = norm_b * (mid_point_b / norm(mid_point_b))
        point_c = norm_c * (mid_point_c / norm(mid_point_c))
        point_d = norm_d * (mid_point_d / norm(mid_point_d))
        #
        # Weight the control points
        weight_int = norm(mid_point_a) / r_int
        weight_ext = norm(mid_point_c) / r_ext

        print(f"WEIGHT INT: {r_int}")
        print(f"WEIGHT EXT: {r_ext}")
        # NOTE: the computations were done in 3D, but the weighted control
        # points are stored in 4D.
        point_a = (weight_int * point_a).tolist() + [weight_int*corner_1[2], weight_int]
        point_b = (weight_int * point_b).tolist() + [weight_int*corner_2[2], weight_int]
        point_c = (weight_ext * point_c).tolist() + [weight_ext*corner_5[2], weight_ext]
        point_d = (weight_ext * point_d).tolist() + [weight_ext*corner_6[2], weight_ext]

        if flat_face == FlatFace.INT:
            #
            # The face 1-2-3-4 is flat
            point_a = (mid_point_a).tolist() + [corner_1[2], 1.]
            point_b = (mid_point_b).tolist() + [corner_2[2], 1.]
            #
        elif flat_face == FlatFace.EXT:
            #
            # The face 5-6-7-8 is flat
            point_c = (mid_point_c).tolist() + [corner_5[2], 1.]
            point_d = (mid_point_d).tolist() + [corner_6[2], 1.]

        control_points = asarray([
            [
                [
                    corner_1.tolist(),
                    point_a,
                    corner_4.tolist()
                ],
                [
                    corner_5.tolist(),
                    point_c,
                    corner_8.tolist()
                ],
            ],
            [
                [
                    corner_2.tolist(),
                    point_b,
                    corner_3.tolist()
                ],
                [
                    corner_6.tolist(),
                    point_d,
                    corner_7.tolist()
                ],
            ]
        ])

        return control_points
