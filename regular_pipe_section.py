from general_pipe_section import GeneralPipeSection
from numpy import sin, cos, asarray
from enums import FlatFace

class RegularPipeSection(GeneralPipeSection):

    def __init__(self,
        planar_rotation=None,
        translation=None,
        *,
        inner_radius,
        outer_radius,
        height,
        angle
    ):

        corner_2 = asarray([inner_radius*cos(angle), inner_radius*sin(angle), 0., 1.])
        corner_3 = asarray([inner_radius,                                  0, 0., 1.])

        corner_6 = asarray([outer_radius*cos(angle), outer_radius*sin(angle), 0., 1.])
        corner_7 = asarray([outer_radius,                                  0, 0., 1.])

        vertical_offset = asarray([0, 0, height, 0])
        corner_1 = corner_2 + vertical_offset
        corner_4 = corner_3 + vertical_offset
        corner_5 = corner_6 + vertical_offset
        corner_8 = corner_7 + vertical_offset

        super().__init__(
            planar_rotation=planar_rotation,
            translation=translation,
            corners=[
                corner_1,
                corner_2,
                corner_3,
                corner_4,
                corner_5,
                corner_6,
                corner_7,
                corner_8,
            ], 
            flat_face=FlatFace.NONE
        )
        
        return