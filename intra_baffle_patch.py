from general_pipe_section import GeneralPipeSection
from enums import FlatFace
from numpy import sqrt, asarray, cos, sin, dot
from numpy.linalg import norm

class IntraBafflePatch(GeneralPipeSection):

    def __init__(self,
        translation=None,
        *,
        intra_baffle_angle,
        baffle_length,
        baffle_depth,
        baffle_width,
        baffle_number,
        vessel_radius
    ):
    
        #
        # Fancy trig
        y = baffle_width/2.
        x = sqrt(vessel_radius**2. - (baffle_width/2.)**2.)
        corner_6 = self.rotate(asarray([x, -y, 0, 1.]), intra_baffle_angle)
        corner_7 = asarray([x,  y, 0, 1.])

        #
        # Extend straight inward in the r-direction
        corner_2 = self.rotate(asarray([x-baffle_depth, -y, 0, 1.]), intra_baffle_angle)
        corner_3 = asarray([x-baffle_depth, y, 0, 1.])
        self.r_int = norm(corner_2[:2])

        #
        # Extrude in the z-direction
        extrusion = asarray([0, 0, baffle_length, 0])
        corner_1 = corner_2 + extrusion
        corner_4 = corner_3 + extrusion
        corner_5 = corner_6 + extrusion
        corner_8 = corner_7 + extrusion

        super().__init__(
            planar_rotation=intra_baffle_angle*baffle_number,
            translation=translation,
            corners=[
                corner_1,
                corner_2,
                corner_3,
                corner_4,
                corner_5,
                corner_6,
                corner_7,
                corner_8
            ],
            flat_face=FlatFace.NONE
        )
        return

    def rotate(self, point, angle):
        Rz = asarray([
            [cos(angle), -sin(angle), 0],
            [sin(angle),  cos(angle), 0],
            [         0,           0, 1]
        ])
        point[:3] = dot(Rz, point[:3])
        return point