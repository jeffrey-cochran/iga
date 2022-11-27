from volumetric_union import VolumetricUnion
from shaft_patches import ShaftPatches
from numpy import pi, asarray
from typing import List
from enums import View

class InternalVesselVolume(VolumetricUnion):

    def __init__(self,
        view=View.NONE,
        *,
        num_impellers:int,
        num_blades_list:List[int],
        blade_length_list:List[float],
        blade_width_list:List[float],
        blade_depth_list:List[float],
        gaps_list:List[float],
        shaft_radius:float,
        bounding_radius:float,
    ):

        #
        # Check compatibility
        assert(len(num_blades_list) == num_impellers)
        assert(len(blade_length_list) == num_impellers)
        assert(len(blade_width_list) == num_impellers)
        assert(len(blade_depth_list) == num_impellers)
        assert(len(gaps_list) == num_impellers)
        assert(bounding_radius > shaft_radius)

        volumetric_patches = []

        accumulated_z_translation = 0
        for i in range(num_impellers-1,-1,-1):
            impeller_translation=asarray([0.,0.,accumulated_z_translation])
            if view==View.STAGGER:
                impeller_translation+=asarray([2*bounding_radius,0.,0.])
            ip = ImpellerPatches(
                translation=impeller_translation,
                num_blades=num_blades_list[i],
                blade_depth=blade_depth_list[i],
                blade_length=blade_length_list[i],
                blade_width=blade_width_list[i],
                shaft_radius=shaft_radius,
                bounding_radius=bounding_radius
            )
            #
            volumetric_patches += ip.volumetric_patches
            accumulated_z_translation += blade_length_list[i]

            sp = ShaftPatches(
                translation=asarray([0,0,accumulated_z_translation]),
                shaft_length=gaps_list[i],
                shaft_radius=shaft_radius,
                bounding_radius=bounding_radius
            )
            #
            volumetric_patches += sp.volumetric_patches
            accumulated_z_translation += gaps_list[i]

        super().__init__(
            volumetric_patches=volumetric_patches,
            offset_v=0,
            offset_vp=0,
            offset_trim=0
        )

        return

