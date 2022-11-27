from intra_blade_patch import IntraBladePatch
from regular_pipe_section import RegularPipeSection
from volumetric_union import VolumetricUnion
from numpy import pi, sin

class ImpellerPatches(VolumetricUnion):

    def __init__(self,
        translation=None,
        *,
        num_blades,
        blade_width,
        blade_length,
        blade_depth,
        shaft_radius,
        bounding_radius
    ):

        intra_blade_angle= 2*pi/num_blades

        assert(sin(intra_blade_angle) >= blade_width / (2*shaft_radius))
        assert(1 > blade_width / (2*shaft_radius))

        blade_patches = [
            IntraBladePatch(
                translation=translation,
                intra_blade_angle=intra_blade_angle,
                blade_length=blade_length,
                blade_depth=blade_depth,
                blade_width=blade_width,
                blade_number=i,
                shaft_radius=shaft_radius
            ) for i in range(num_blades)
        ]

        self.pseudo_radius = blade_patches[0].r_ext
        
        assert(bounding_radius > self.pseudo_radius)
        bounding_cylinder_patches = [
            RegularPipeSection(
                translation=translation,
                planar_rotation=pi*(i/2.),
                inner_radius=self.pseudo_radius,
                outer_radius=bounding_radius,
                height=blade_length,
                angle=pi/2.
            ) for i in range(4)
        ]
            
        combined_patches =  blade_patches + bounding_cylinder_patches

        super().__init__(volumetric_patches=combined_patches)
