from patch import Patch
from numpy import typing, asarray

class TrimCurve(Patch):

    def __init__(
        self, *,
        knots,
        point_1,
        point_2,
        degree,
        offset
    ):

        self.u0 = knots[0]
        self.u1 = knots[-1]
        self.knots = knots
        self.point_1 = point_1
        self.point_2 = point_2
        self.degree = degree
        self.offset = offset
        
        return

    def dumps(self):  
        knot_string = [str(x) for x in self.knots]
        knot_string = " ".join(knot_string)          
        return f"""# Begin trim curve
vp {self.point_1[0]} {self.point_1[1]}
vp {self.point_2[0]} {self.point_2[1]}
cstype bspline
deg {self.degree}
curv2 {self.offset+1} {self.offset+2}
parm u {knot_string}
end

"""