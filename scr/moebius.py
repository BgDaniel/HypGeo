from complex_plane import *

class Moeb:
    """A class representing Moebius transformations (of the connected unit component) and their operation on the upper half plane."""
    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def d(self):
        return self._d

    def __init__(self, a, b, c, d):
        """
        Parameters
        ----------
        a : float or (real) ComplexNumer
            coefficient c_11 in matrix representation of moebius transformation via SL(2, R)
            given by ( (c_11, c_12), (c_21, c_22))
        b : float or (real) ComplexNumerr
            coefficient c_12 in matrix representation
        c : float or (real) ComplexNumerr
            coefficient c_21 in matrix representation
        a : float or ComplexNumer
            coefficient c_22 in matrix representation
        """
        if type(a) is float:
            a = ComplexNumber(a, .0)
        if type(b) is float:
            b = ComplexNumber(b, .0)
        if type(c) is float:
            c = ComplexNumber(c, .0)
        if type(d) is float:
            d = ComplexNumber(d, .0)
        self._a = a
        self._b = b
        self._c = c
        self._d = d

    def __call__(self, z):
        if type(z) is np.array or type(z) is np.ndarray:
            assert len(z) == 2, "z has wrong dimension!"
            z = ComplexNumber(z[0], z[1])
        elif type(z) is not ComplexNumber:
            raise Exception("Argument z is neither of type array nor type ComplexNumber!") 

        return (z * self._a + self._b) / (z * self._c + self._d)

    def __mul__(self, o):
        """Right multiplicates of instance transformation with o.

        Parameters
        ----------
        o : Moeb
            Moebius transformation multiplied with instance from right

        Returns
        -------
        Moeb
            Right multiplication of instance transformation with o
        """ 
        assert type(o) is Moeb or type(o) is MoebConj, "%s is not of type \"Moeb\"!" % str(o)

        return Moeb(self._a * o.a + self._b * o.c, self._a * o.b + self._b * o.d, 
            self._c * o.a + self._d * o.c, self._c * o.b + self._d * o.d)

    def conj(self, o):
        """Returns the (group) conjugation of instance transformation with o.

        Parameters
        ----------
        o : Moeb
            Moebius transformation used in conjugation of instance.

        Returns
        -------
        Moeb
            (Group) conjugation of instance transformation by multiplying o from the right and its inverse
            from the left
        """ 
        assert type(o) is Moeb, "%s is not of type \"Moeb\"!" % str(o)

        return (o.inverse()) * self * o

    def inv(self):
        """Returns the (group) inverse of instance."""
        return Moeb(self._d, - self._b, - self._c, self._a)

    def plot(self, rectanlge):
        return None

class MoebConj(Moeb):
    """A class representing Moebius transformations (of the second connected component) and their operation on the upper half plane."""
    def __init__(self, a, b, c, d):
        Moeb.__init__(self, a, b, c, d)

    def __call__(self, z):
        if type(z) is np.array or type(z) is np.ndarray or type(z) is list:
            assert len(z) == 2, "z has wrong dimension!"
            z = ComplexNumber(z[0], z[1])
        elif type(z) is not ComplexNumber:
            raise Exception("Argument z is neither of type array nor type ComplexNumber!") 
        return Moeb.__call__(self, z.conj())

class MoebCorr(Moeb):
    def __init__(self, a, b, c, d, rho):
        self._rho = rho
        Moeb.__init__(self, a, b, c, d)

    def __call__(self, z):
        z_x = z.real_part
        z_y = z.imaginary_part
        u = z_x / self._rho + math.sqrt(1.0 - self._rho * self._rho) / self._rho * z_y
        v = z_y
        z = ComplexNumber(u, v)
        w = Moeb.__call__(self, z)
        w_x = w.real_part
        w_y = w.imaginary_part
        u = self._rho * w_x + math.sqrt(1.0 - self._rho * self._rho) * w_y
        v = w_y
        return ComplexNumber(u, v)






