import random
from math import ceil
from finite_field import FiniteFieldElement, FiniteField


class EllipticCurve:
    """
    Class representing an elliptic curve
    in the simplified form.
    Inspired by: https://github.com/j2kun/elliptic-curves-finite-fields
    """
    def __init__(self, a, b, ff):
        """
        Elliptic curve in the simplified form.
        y^2 = x^3 + ax + b
        :param a: Coefficient a
        :param b: Coefficient b
        :param ff: FiniteField of the curve
        """
        if not isinstance(a, FiniteFieldElement) or not isinstance(b, FiniteFieldElement):
            raise TypeError('This elliptic curve only supports finite field elements')
        self.a = a
        self.b = b
        self.finite_field = ff

        if 4*a**3 + 27*b**2 == 0:
            raise ValueError('Bad curve: 4a^3 + 27b^2 mustn\'t be 0')

    def is_point(self, x, y):
        """
        Test whether the point [x,y]
        belongs to the curve.
        :param x: X coordinate
        :param y: Y cooordinate
        :return: bool
        """
        return y ** 2 == x ** 3 + self.a * x + self.b

    def __str__(self):
        return f'y^2 = x^3 + {self.a}x + {self.b}'

    def random_point(self):
        x = 0
        res = 0
        while True:
            x = FiniteFieldElement(random.randint(0, self.finite_field.modulo), self.finite_field.modulo)
            res = x**3 + self.a * x + self.b
            quadr = res ** ((self.finite_field.modulo-1)/2)
            if quadr == 1:
                break

        y_2 = FiniteFieldElement(res, self.finite_field.modulo)
        y = y_2.square_root()
        return ECPoint(x, y, self)

    def order(self):
        while True:
            point = self.random_point()
            m = ceil(self.finite_field.modulo ** float(1/4))

            baby_steps = []
            for j in range(m+1):
                baby_steps.append(j*point)

            q = (2*m + 1)*point
            r = (self.finite_field.modulo+1)*point

            matches = []
            t = ceil((2 * self.finite_field.modulo**float(1/2))/(2*m+1))

            for i in range(1, t+1):
                new_point = r + i*q
                for j in range(len(baby_steps)):
                    if (new_point.x, new_point.y) == (baby_steps[j].x, baby_steps[j].y):
                        matches.append(self.finite_field.modulo + 1 + i*(2*m+1)-j)
            if len(matches) == 1:
                return matches[0]
            print(f'Unsuccessful bsgs: {len(matches)} matches')


class ECPoint:
    """
    Class representing a point on
    an elliptic curve
    """

    def __init__(self, x, y, curve):
        self.curve = curve
        self.x = x
        self.y = y

        if not curve.is_point(x, y):
            raise ValueError(f'The point {self} is not on the curve {curve}!')

    def __str__(self):
        return f'[{self.x}, {self.y}]'

    def __neg__(self):
        """
        Negation of the point.
        -P = (x,-y)
        :return: ECPoint
        """
        return ECPoint(self.x, - self.y, self.curve)

    def __add__(self, other):
        """
        Add operation of two EC points
        as defined in the handout.
        :param other: Other EC point
        :return: ECPoint
        """

        if self.curve != other.curve:
            raise ValueError('Can\'t add points on different curves!')

        # P + 0 = P
        if isinstance(other, ECPointAtInfinity):
            return self

        # 0 + P = P
        if isinstance(self, ECPointAtInfinity):
            return other

        # x1 = x2 AND y1 = -y2 => P + Q = 0
        if (self.x, self.y) == (other.x, -other.y):
            return ECPointAtInfinity(self.curve)

        # Prepare lambda parameter
        l = 0

        # P = Q
        if (self.x, self.y) == (other.x, other.y):
            try:
                l = (3 * (self.x ** 2) + self.curve.a) / (2 * self.y)
            except ZeroDivisionError:
                return ECPointAtInfinity(self.curve)
        # P != Q
        else:
            try:
                l = (other.y - self.y) / (other.x - self.x)
            # x1 == x2
            except ZeroDivisionError:
                return ECPointAtInfinity(self.curve)

        r1 = l**2 - self.x - other.x
        r2 = l*(self.x - other.x) - self.y
        return ECPoint(r1, r2, self.curve)

    def __sub__(self, other):
        """
        Substraction of two EC points.
        :param other: Other EC point
        :return: ECPoint
        """
        return self + (-other)

    def __mul__(self, n):
        """
        Calculate A = n*B using double-and-add algorithm.
        :param n: Factor of multiplication
        :return: ECPoint
        """
        if not isinstance(n, int):
            raise TypeError(f'Can\'t multiply an ECPoint by {type(n)}!')
        else:
            if n < 0:
                return -self * -n
            if n == 0:
                return ECPointAtInfinity(self.curve)
            Q = self
            R = self if n & 1 == 1 else ECPointAtInfinity(self.curve)
            i = 2
            while i <= n:
                # double
                Q = Q + Q
                if n & i == i:
                    # add
                    R = Q + R
                i = i << 1
            return R

    def __rmul__(self, n):
        """
        Multiplication with point as the value on the right.
        :param n: Factor
        :return: ECPoint
        """
        return self * n


class ECPointAtInfinity(ECPoint):
    """
    Special case of "zero" point on the elliptic curve.
    """

    def __init__(self, curve):
        self.curve = curve

    def __neg__(self):
        """
        Negated PaI is PaI
        :return: ECPointAtInfinity
        """
        return self

    def __str__(self):
        return 'Point at infinity'

    def __add__(self, Q):
        """
        PaI + Q = Q
        :param Q: Other point
        :return: Point Q
        """
        if self.curve != Q.curve:
            raise ValueError('Can\'t add points on different curves!')
        return Q

    def __mul__(self, n):
        """
        Multiplied PaI is PaI
        :param n: Multiplication factor
        :return: ECPointAtInfinity
        """
        if not isinstance(n, int):
            raise TypeError(f'Can\'t multiply a point by {type(n)}!')
        return self

    def __eq__(self, other):
        """
        Test whether other point is also PaI
        :param other: Other point
        :return: bool
        """
        return type(other) is ECPointAtInfinity

