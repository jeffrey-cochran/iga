from patch import Patch
from volumetric_patch import VolumetricPatch
from typing import List, AnyStr
from numpy import asarray
import tempfile

class VolumetricUnion(Patch):

    def __init__(
        self,
        volumetric_patches: List[VolumetricPatch]=[],
        offset_v=0,
        offset_vp=0,
        offset_trim=0
    ):
        self.volumetric_patches=volumetric_patches

        self.current_offset_v = offset_v
        self.current_offset_vp = offset_vp
        self.current_offset_trim = offset_trim

        if len(volumetric_patches) > 0:
            self.compute_offsets()

        return

    def __add__(self, other_volumetric_union):
        return VolumetricUnion(
            volumetric_patches = self.volumetric_patches + other_volumetric_union.volumetric_patches,
            offset_v=0,
            offset_vp=0,
            offset_trim=0
        )

    def compute_offsets(self):
        for volumetric_patch in self.volumetric_patches:
            self.update_offsets(volumetric_patch)
        return

    def update_offsets(self, volumetric_patch:VolumetricPatch):
        volumetric_patch.set_offsets(
            offset_v=self.current_offset_v, 
            offset_vp=self.current_offset_vp, 
            offset_trim=self.current_offset_trim
        )
        self.current_offset_v += volumetric_patch.num_v
        self.current_offset_vp += volumetric_patch.num_vp
        self.current_offset_trim += volumetric_patch.num_trimming_curves
        return

    def append_volumetric_patch(self, volumetric_patch:VolumetricPatch):
        self.update_offsets(volumetric_patch)
        self.volumetric_patches.append(volumetric_patch)
        return

    def dumps(self):
        volume_str = ""
        for volumetric_patch in self.volumetric_patches:
            volume_str += volumetric_patch.dumps()

        return f"""# Begin volume union dump
{volume_str}
"""