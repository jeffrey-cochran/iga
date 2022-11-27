from patch import Patch
from trim_curve import TrimCurve
from surface_patch import SurfacePatch
from typing import List, AnyStr
from numpy import asarray
import tempfile

class SurfaceUnion(Patch):

    def __init__(
        self,
        surface_patches: List[SurfacePatch]=[],
        offset_v=0,
        offset_vp=0,
        offset_trim=0, 
        *,
        object_name: AnyStr,
        layer_name: AnyStr
    ):
        self.surface_patches=surface_patches
        self.object_name=object_name
        self.layer_name=layer_name

        self.init_offset_v = offset_v
        self.init_offset_vp = offset_vp
        self.init_offset_trim = offset_trim

        self.current_offset_v = offset_v
        self.current_offset_vp = offset_vp
        self.current_offset_trim = offset_trim

        if len(surface_patches) > 0:
            self.compute_offsets()

        return

    def compute_offsets(self):
        for surface_patch in self.surface_patches:
            self.update_offsets(surface_patch)
        return

    def update_offsets(self, surface_patch:SurfacePatch):
        surface_patch.set_offsets(
            offset_v=self.current_offset_v, 
            offset_vp=self.current_offset_vp, 
            offset_trim=self.current_offset_trim
        )
        self.current_offset_v += surface_patch.num_v
        self.current_offset_vp += surface_patch.num_vp
        self.current_offset_trim += surface_patch.num_trim_curves
        return

    def object_counts(self):
        return (
            self.current_offset_v - self.init_offset_v,
            self.current_offset_vp - self.init_offset_vp,
            self.current_offset_trim - self.init_offset_trim
        )

    def append_surface_patch(self, surface_patch):
        self.update_offsets(surface_patch)
        self.surface_patches.append(surface_patch)
        return

    def dumps(self, include_trim_loop=True):
        surfaces_str = ""
        for surface in self.surface_patches:
            surfaces_str += surface.dumps()

        return f"""# Begin surface union dump
g {self.layer_name}
o {self.object_name}

{surfaces_str}
"""