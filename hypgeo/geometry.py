import enum
from hypgeo.complex_plane import *
import math

class LineType(enum.Enum):
    CIRCLE = 1
    VERTICAL = 2

class GeodesicLine:
    def line_trough(z0, z1):
        assert z0 != z1, "z0 and z1 have to be different in order to determine a geodesic line in upper half plane!"
        assert type(z0) is ComplexNumber, "z0 is not of type \"ComplexNumber\"!"
        assert type(z1) is ComplexNumber, "z1 is not of type \"ComplexNumber\"!"
        assert z0 != z1, "A line in HPlus can only be determined by two different complex numbers!"
        assert z0.im > 0, "z0 is not contained in HPlus!"
        assert z1.im > 0, "z1 is not contained in HPlus!"

        if math.isclose(z0.re, z1.re, abs_tol=1e-09):
            return Vertical(z0.re)
        else:
            x0 = self._z0.re
            y0 = self._z0.im
            x1 = self._z1.re
            y1 = self._z1.im
            center = 1.0 / 2.0 * (y1 ** 2 - y0 ** 2 - x0 ** 2 + x1 ** 2) / (x1 - x0)
            radius = 1.0 / ( 2.0 * (abs(x1 - x0))) * math.sqrt(((x1 - x0) ** 2 + (y1 - y0) ** 2) * ((x1 - x0) ** 2 + (y1 + y0) ** 2))
            return HalfCircle(radius, center)

    def __init__(self):
        pass

class Vertical(GeodesicLine):
    @property
    def Absc(self):
        return self._absc

    def __init__(self, absc):
        self._absc = absc

    def __eq__(self, o):
        if type(o) != Vertical:
            return False
        else:
            if math.isclose(self._absc, o.Absc, abs_tol=1e-09):
                return True
            else: 
                return False

    def rnd_v(min_absc, max_absc, samples):
        rnd_v = []
        Absc = np.random.uniform(min_absc, max_absc, samples)

        for absc in Absc:       
            rnd_v.append(Vertical(absc))

        return rnd_v
    
class HalfCircle(GeodesicLine):
    @property
    def Radius(self):       
        return self._radius

    @property
    def Center(self):
        return self._center

    @property
    def Level(self):
        return self._level

    def __init__(self, radius, center):
        self._radius = radius
        self._center = center

    def __eq__(self, o):
        if type(o) != HalfCircle:
            return False
        else:
            if math.isclose(self._radius, o.Radius, abs_tol=1e-09) and math.isclose(self._center, o.Center, abs_tol=1e-09):
                return True
            else: 
                return False

    def rnd_c(max_r, min_c, max_c, samples):
        rnd_c = []
        R, C = np.random.uniform(.0, max_r, samples), np.random.uniform(min_c, max_c, samples)

        for i in range(0, samples):
            rnd_c.append(HalfCircle(R[i], C[i]))


            HalfCircle
            z_0, z_1 = ComplexNumber(C[i] + R[i] * math.sin(math.pi / 4.0), R[i] * math.cos(math.pi / 4.0)), \
                ComplexNumber(C[i] - R[i] * math.sin(math.pi / 4.0), R[i] * math.cos(math.pi / 4.0))        
            rnd_c.append(Line(z_0, z_1))

        return rnd_c

unit_circle = HalfCircle(1.0, .0)

class Position(enum.Enum):
    IN = 1
    OUT = 2
    ON = 3

def opposite(position):
    if position == Position.IN:
        return Position.OUT
    elif position == Position.OUT:
        return Position.IN
    else:
        return Position.ON

class HalfSpace:
    def __init__(self, line):
        self._line = line
        self._refl = refl(line)
        
    def position(self, x):
        if type(x) is np.array or type(x) is np.ndarray or type(x) is list:
            assert len(x) == 2, "z has wrong dimension!"
        elif type(x) is ComplexNumber:
            x = x.ToVector()
        else:
            raise Exception("Argument z is neither of type array nor type ComplexNumber!") 
        level = self._line.Level(x[0], x[1])

        if level < .0:
            return Position.IN
        elif level == .0:
            return Position.ON
        else:
            return Position.OUT

    def refl(self, x):
        return self._refl(x) 
