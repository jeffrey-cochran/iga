from internal_vessel_volume import InternalVesselVolume
from baffle_patches import BafflePatches
from numpy import pi, asarray
from enums import View

num_blades = 6
intra_blade_angle= 2*pi/num_blades
blade_length=1.
blade_width=blade_length*.05
blade_depth=blade_length*.5
shaft_radius= 0.1
bounding_radius = 1.5*(shaft_radius+blade_depth)

ivv = InternalVesselVolume(
    translation=asarray([0.,0,0]),
    view=View.NONE,
    num_impellers=2,
    num_blades_list=[6,6],
    blade_depth_list=[blade_depth,blade_depth],
    blade_width_list=[blade_width,blade_width],
    blade_length_list=[blade_length,blade_length],
    gaps_list=[3.,3.],
    shaft_radius=shaft_radius,
    bounding_radius=bounding_radius
)

bp = BafflePatches(
    translation=None,
    num_baffles=4,
    baffle_width=0.05,
    baffle_length=8,
    baffle_depth=.5,
    vessel_radius=8/3.,
    bounding_radius=bounding_radius
)

combo = ivv + bp
with open("pipe_section.obj", "w") as f:
    f.write(combo.dumps())

