from regular_pipe_section import RegularPipeSection
from intra_baffle_patch import IntraBafflePatch
from volumetric_union import VolumetricUnion
from numpy import pi, sin

class BafflePatches(VolumetricUnion):

    def __init__(self,
        translation=None,
        *,
        num_baffles,
        baffle_width,
        baffle_length,
        baffle_depth,
        vessel_radius,
        bounding_radius
    ):

        intra_baffle_angle= 2*pi/num_baffles

        assert(sin(intra_baffle_angle) >= baffle_width / (2*vessel_radius))
        assert(1 > baffle_width / (2*vessel_radius))

        baffle_patches = [
            IntraBafflePatch(
                translation=translation,
                intra_baffle_angle=intra_baffle_angle,
                baffle_length=baffle_length,
                baffle_depth=baffle_depth,
                baffle_width=baffle_width,
                baffle_number=i,
                vessel_radius=vessel_radius
            ) for i in range(num_baffles)
        ]

        self.pseudo_radius = baffle_patches[0].r_int
        
        assert(bounding_radius < self.pseudo_radius)
        bounding_cylinder_patches = [
            RegularPipeSection(
                translation=translation,
                planar_rotation=pi*(i/2.),
                inner_radius=bounding_radius,
                outer_radius=self.pseudo_radius,
                height=baffle_length,
                angle=pi/2.
            ) for i in range(4)
        ]
            
        combined_patches =  baffle_patches + bounding_cylinder_patches

        super().__init__(volumetric_patches=combined_patches)
