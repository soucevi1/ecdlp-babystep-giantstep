
# Module for elliptic curve operations.
# Author: Vit Soucek

from math import inf, sqrt
from finite_field import FiniteFieldElement


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
        if isinstance(a, int):
            a = ff.get_element(a)
        if isinstance(b, int):
            b = ff.get_element(b)
        self.a = a
        self.b = b
        self.finite_field = ff

        if 4 * a ** 3 + 27 * b ** 2 == 0:
            raise ValueError('Bad curve: 4a^3 + 27b^2 mustn\'t be 0')

    def is_point(self, x, y):
        """
        Test whether the point [x,y]
        belongs to the curve.
        :param x: X coordinate
        :param y: Y coordinate
        :return: bool
        """
        if isinstance(x, int):
            x = self.finite_field.get_element(x)
        if isinstance(y, int):
            y = self.finite_field.get_element(y)
        return y ** 2 == x ** 3 + self.a * x + self.b

    def __str__(self):
        return f'y^2 = x^3 + {self.a}x + {self.b}'

    def __eq__(self, other):
        """
        Overloaded == operator
        :param other: Other elliptic curve
        :return: Bool
        """
        return (self.a, self.b) == (other.a, other.b)

    def order_approx(self):
        """
        Approximation according to Hasse theorem.
        :return: Upper bound of the order.
        """
        return self.finite_field.modulo + 1 + 2 * sqrt(self.finite_field.modulo)


class ECPoint:
    """
    Class representing a point on
    an elliptic curve
    """

    def __init__(self, x, y, curve):
        self.curve = curve
        self.x = FiniteFieldElement(x, self.curve.finite_field.modulo)
        self.y = FiniteFieldElement(y, self.curve.finite_field.modulo)

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
        as defined in the handouts.
        :param other: Other EC point
        :return: ECPoint
        """

        if self.curve != other.curve:
            raise ValueError('Cannot add points on different curves!')

        # P + 0 = P
        if isinstance(other, ECPointAtInfinity):
            return self

        # 0 + P = P
        if isinstance(self, ECPointAtInfinity):
            return other

        # x1 = x2 AND y1 = -y2 => P + Q = 0
        if (self.x, self.y) == (other.x, -other.y):
            return ECPointAtInfinity(self.curve)

        # P = Q
        if (self.x, self.y) == (other.x, other.y):
            try:
                lam = (3 * (self.x ** 2) + self.curve.a) / (2 * self.y)
            except ZeroDivisionError:
                return ECPointAtInfinity(self.curve)
        # P != Q
        else:
            try:
                lam = (other.y - self.y) / (other.x - self.x)
            # x1 == x2
            except ZeroDivisionError:
                return ECPointAtInfinity(self.curve)

        r1 = lam ** 2 - self.x - other.x
        r2 = lam * (self.x - r1) - self.y
        return ECPoint(r1, r2, self.curve)

    def __sub__(self, other):
        """
        Difference of two EC points.
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
            raise TypeError(f'Cannot multiply an ECPoint by {type(n)}!')
        else:
            if n < 0:
                return -self * -n
            if n == 0:
                return ECPointAtInfinity(self.curve)
            point_q = self
            point_r = self if n & 1 == 1 else ECPointAtInfinity(self.curve)
            i = 2
            while i <= n:
                # double
                point_q = point_q + point_q
                if n & i == i:
                    # add
                    point_r = point_q + point_r
                i = i << 1
            return point_r

    def __rmul__(self, n):
        """
        Multiplication with point as the value on the right.
        :param n: Factor
        :return: ECPoint
        """
        return self * n

    def __eq__(self, other):
        """
        Overloaded == operator
        :param other: Other ECPoint
        :return: True if self == other, False otherwise
        """
        return (self.x, self.y, self.curve) == (other.x, other.y, other.curve)

    def order_approx(self):
        """
        The upper bound of the EC point order
        is the order of the EC itself.
        :return: Approx. of the order of the EC
        """
        return self.curve.order_approx()

    def __lt__(self, other):
        """
        Overloaded less-than operator
        for sorting purposes.
        :param other: Other ECPoint
        :return: Bool
        """
        if self.x < other.x:
            return True
        if self.x > other.x:
            return False
        return self.y < other.y


class ECPointAtInfinity(ECPoint):
    """
    Special case of "zero" point on the elliptic curve.
    """

    def __init__(self, curve):
        self.x = inf
        self.y = inf
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
            raise TypeError(f'Cannot multiply a point by {type(n)}!')
        return self

    def __eq__(self, other):
        """
        Test whether other point is also PaI
        :param other: Other point
        :return: bool
        """
        return type(other) is ECPointAtInfinity
