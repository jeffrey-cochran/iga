from intra_blade_patch import IntraBladePatch
from regular_pipe_section import RegularPipeSection
from volumetric_union import VolumetricUnion
from numpy import pi, sin

class ShaftPatches(VolumetricUnion):

    def __init__(self,
        translation=None,
        *,
        shaft_length,
        shaft_radius,
        bounding_radius,
    ):
        
        assert(bounding_radius > shaft_radius)

        bounding_cylinder_patches = [
            RegularPipeSection(
                translation=translation,
                planar_rotation=pi*(i/2.),
                inner_radius=shaft_radius,
                outer_radius=bounding_radius,
                height=shaft_length,
                angle=pi/2.
            ) for i in range(4)
        ]

        super().__init__(volumetric_patches=bounding_cylinder_patches)
