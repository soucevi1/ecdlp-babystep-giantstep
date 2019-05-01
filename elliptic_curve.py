import random
from math import ceil, inf, gcd, sqrt
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
        if isinstance(a, int):
            a = ff.get_element(a)
        if isinstance(b, int):
            b = ff.get_element(b)
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
        if isinstance(x, int):
            x = self.finite_field.get_element(x)
        if isinstance(y, int):
            y = self.finite_field.get_element(y)
        return y ** 2 == x ** 3 + self.a * x + self.b

    def __str__(self):
        return f'y^2 = x^3 + {self.a}x + {self.b}'

    def random_point(self):
        x = 0
        res = 0
        p = int((self.finite_field.modulo-1)/2)
        while True:
            x = FiniteFieldElement(random.randint(0, self.finite_field.modulo), self.finite_field.modulo)
            res = x**3 + self.a * x + self.b
            quadr = res ** p
            if quadr == 1:
                break

        y_2 = FiniteFieldElement(res, self.finite_field.modulo)
        y = y_2.square_root()
        return ECPoint(x, y, self)

    def order_approx(self):
        """
        Approximation according to Hasse theorem.
        :return: Upper bound of the order.
        """
        return self.finite_field.modulo + 1 + 2*sqrt(self.finite_field.modulo)


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
        r2 = l*(self.x - r1) - self.y
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

    def __eq__(self, other):
        """
        Overloaded == operator
        :param other: Other ECPoint
        :return: True if self == other, False otherwise
        """
        return (self.x, self.y) == (other.x, other.y)

    def order(self):
        """
        Get an order of the EC point
        using babystep-giantstep algorithm.
        Inspired by: https://en.wikipedia.org/wiki/Counting_points_on_elliptic_curves#Baby-step_giant-step
        :return: The order of the point
        """
        # Normal start of BSGS to find order of the whole elliptic curve
        m = ceil(self.curve.finite_field.modulo ** float(1 / 4))
        P = self
        baby_steps = []

        # Prepare the baby steps
        for j in range(m + 1):
            baby_steps.append(j * P)

        Q = (self.curve.finite_field.modulo + 1) * P

        k = 0
        j = 0
        while True:
            # Calculate one giant step
            new_point = Q + k * (2 * m * P)
            br = False

            # Compare the GS with all BS to find a match
            for i in range(len(baby_steps)):
                # Compare with j*P and also -j*P
                if new_point == baby_steps[i]:
                    j = -i
                    br = True
                    break
                if new_point == -baby_steps[i]:
                    j = i
                    br = True
            if br:
                break
            k += 1

        # Calculate M such that M*P = 0
        # => M is a multiple of the order of P
        M = self.curve.finite_field.modulo + 1 + 2 * m * k + j

        # Factor M to prime factors
        factors = factor(M)

        # Divide M by its factors
        # until the smallest possible x
        # such as x*P = 0 is found
        for f in factors:
            x = int(M / f)
            xp = x * P
            if isinstance(xp, ECPointAtInfinity):
                M = x

        return M

    def order_approx(self):
        """
        The upper bound of the EC point order
        is the order of the EC itself.
        :return: Approx. of the order of the EC
        """
        return self.curve.order_approx()

    def __lt__(self, other):
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
            raise TypeError(f'Can\'t multiply a point by {type(n)}!')
        return self

    def __eq__(self, other):
        """
        Test whether other point is also PaI
        :param other: Other point
        :return: bool
        """
        return type(other) is ECPointAtInfinity


def factor(n):
    """
    Factor number to prime factors.
    Source: https://stackoverflow.com/a/22808285/6136143
    :param n: Number to factor
    :return: List of factors.
    """
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors


def lcm(a, b):
    """
    Find the least common multiple of a and b.
    Source: https://stackoverflow.com/a/51716959/6136143
    :param a: First number
    :param b: Second number
    :return: Least common multiple of a and b
    """
    return abs(a*b) // gcd(a, b)


def binary_search(element_list, element):
    """
    Source: https://stackoverflow.com/a/4161779/6136143
    :param element_list: List of elements (BabyStepPoints)
    :param element: Element to look for (ECPoint)
    :return:
    """
    hi = len(element_list)
    lo = 0
    while lo < hi:
        mid = (lo+hi)//2
        midval = element_list[mid].point
        if midval < element:
            lo = mid+1
        elif midval > element:
            hi = mid
        else:
            return element_list[mid].index
    return -1
