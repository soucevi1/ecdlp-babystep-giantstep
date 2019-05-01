from math import inf


class FiniteFieldElement:
    """
    Class representing one element of a finite field.
    Every element has its value, modulo and
    special operations on its field.
    """

    def __init__(self, value, modulo):
        """
        Initialize the finite field element with its
        value and modulo of the finite field.
        :param value: Value of the element
        :param modulo: Modulo of the finite field
        """
        self.modulo = modulo
        if isinstance(value, FiniteFieldElement):
            self.value = value.value % modulo
        elif isinstance(value, int):
            self.value = value % modulo
        else:
            raise TypeError(f'Type {type(value)} cannot be used in FFE constructor')

    def __str__(self):
        """
        Overloaded string operator.
        :return:
        """
        return str(self.value)

    def __add__(self, other):
        """
        Overloaded + operator.
        :param other: Othre FiniteFieldElement
        :return: Sum of the elements
        """
        if isinstance(other, int):
            return FiniteFieldElement(self.value + other, self.modulo)
        elif isinstance(other, FiniteFieldElement):
            return FiniteFieldElement(self.value + other.value, self.modulo)
        else:
            raise TypeError(f'Adding {type(other)} to FiniteFieldElement not implemented')

    def __radd__(self, other):
        """
        Overloaded + operator (rvalue).
        :param other: Othre FiniteFieldElement
        :return: Sum of the elements
        """
        return self + other

    def __sub__(self, other):
        """
        Overloaded - operator.
        :param other: Other FiniteFieldElement
        :return: Difference of the elements
        """
        if isinstance(other, int):
            return FiniteFieldElement(self.value - other, self.modulo)
        elif isinstance(other, FiniteFieldElement):
            return FiniteFieldElement(self.value - other.value, self.modulo)
        else:
            raise TypeError(f'Substracting {type(other)} from FiniteFieldElement not implemented')

    def __rsub__(self, other):
        """
        Overloaded - operator (rvalue).
        :param other: Other FiniteFieldElement
        :return: Difference of the elements
        """
        if isinstance(other, int):
            return FiniteFieldElement(other - self.value, self.modulo)
        elif isinstance(other, FiniteFieldElement):
            return FiniteFieldElement(other.value - self.value, self.modulo)
        else:
            raise TypeError(f'Substracting {type(other)} from FiniteFieldElement not implemented')

    def __mul__(self, other):
        """
        Overloaded * operator.
        :param other: Other FiniteFieldElement
        :return: Multiple of the elements.
        """
        if isinstance(other, int):
            return FiniteFieldElement(self.value * other, self.modulo)
        elif isinstance(other, FiniteFieldElement):
            return FiniteFieldElement(self.value * other.value, self.modulo)
        else:
            raise TypeError(f'Multiplying {type(other)} (value {other}) and FiniteFieldElement not implemented')

    def __rmul__(self, other):
        """
        Overloaded * operator (rvalue).
        :param other: Other FiniteFieldElement
        :return: Multiple of the elements.
        """
        return self * other

    def __neg__(self):
        """
        Negate the element
        :return: Negated element
        """
        return FiniteFieldElement(- self.value, self.modulo)

    def inverse(self):
        """
        Extended Euclidean Algorithm for inverse element.
        Source:
        https://stackoverflow.com/a/17957069/6136143
        :return: Multiplicative inverse of the element.
        """
        if self.value == 1:
            return self.value
        if self.modulo % self.value == 0:
            raise ValueError(f'Failed to invert {self.value}! {self.modulo} is a multiple of {self.value}.')
        a, b = self.modulo, self.value
        x, prevx = 0, 1
        y, prevy = 1, 0
        while b != 0:
            q = a // b
            a, b = b, a % b
            x, prevx = prevx - q * x, x
            y, prevy = prevy - q * y, y
        if a != 1:
            raise ValueError(f'Failed to invert {self.value}! GCD of {self.modulo} and {self.value} is {a}, not 1!')
        return FiniteFieldElement(prevy, self.modulo)

    def __truediv__(self, other):
        """
        Overloaded division operator.
        :param other: Other element or int.
        :return: Division modulo.
        """
        if isinstance(other, int):
            offe = FiniteFieldElement(other, self.modulo)
            iother = offe.inverse()
            return self * iother
        elif isinstance(other, FiniteFieldElement):
            iother = other.inverse()
            return self * iother
        else:
            raise TypeError(f'Division of {type(other)} and FiniteFieldElement not implemented.')

    def __rtruediv__(self, other):
        """
        Overloaded division operator (rvalue).
        :param other: Other element or int.
        :return: Division modulo.
        """
        if isinstance(other, int):
            offe = FiniteFieldElement(other, self.modulo)
            return offe / self
        elif isinstance(other, FiniteFieldElement):
            return other / self
        else:
            raise TypeError(f'Division of {type(other)} and FiniteFieldElement not implemented.')

    def __pow__(self, power):
        """
        Overloaded power operator.
        :param power: Power
        :return: Power of element.
        """
        if isinstance(power, int):
            return FiniteFieldElement(self.value ** power, self.modulo)
        elif isinstance(power, FiniteFieldElement):
            return FiniteFieldElement(self.value ** power.value, self.modulo)
        else:
            raise TypeError(f'Unsupported power type: {type(power)}, value {power}')

    def __rpow__(self, other):
        """
        Overloaded power operator.
        :param other: Power
        :return: Power of element.
        """
        if isinstance(other, int):
            return FiniteFieldElement(other ** self.value, self.modulo)
        elif isinstance(other, FiniteFieldElement):
            return FiniteFieldElement(other.value ** self.value, self.modulo)
        else:
            raise TypeError(f'Unsupported power type {type(other)}')

    def __lt__(self, other):
        """
        Overloaded less-than operator.
        :param other: Other Finite Field Element
        :return: Bool
        """
        if isinstance(other, int):
            return self.value < other
        if isinstance(other, FiniteFieldElement):
            return self.value < other.value
        if other is inf:
            return False
        raise TypeError(f'Comparison between FFE and {type(other)}')

    def __gt__(self, other):
        """
        Overloaded greater-than operator.
        :param other: Other Finite Field Element
        :return: Bool
        """
        if isinstance(other, int):
            return self.value > other
        if isinstance(other, FiniteFieldElement):
            return self.value > other.value
        if other is inf:
            return False
        raise TypeError(f'Comarison between FFE and {type(other)}')

    def square_root(self):
        """
        Get squre root of element over a finite field.
        Source: https://www.researchgate.net/publication/
        320402737_An_Algorithm_to_Find_Square_Root_of_Quadratic_Residues_over_Finite_Fields_using_Primitive_Elements
        :return: FiniteFieldElement
        """
        g = FiniteField(self.modulo).get_primitive_element()
        g0 = g ** 2
        if self.value == g0:
            return g0
        else:
            k = 1
            while (g0 ** k).value != self.value:
                k += 1
            return g ** k

    def __eq__(self, other):
        """
        Overloaded == operator.
        :param other: Other Finite Field Element
        :return: Bool
        """
        if isinstance(other, int):
            return self.value == other
        elif isinstance(other, FiniteFieldElement):
            return self.value == other.value
        elif other is inf:
            return False
        else:
            raise TypeError(f'Cannot compare FF element and {type(other)}')


class FiniteField:
    """
    Class representing an integer
    finite field modulo prime number.
    """

    def __init__(self, modulo):
        """
        Constructor of the finite field with size modulo.
        NOT CHECKING WHETHER modulo IS PRIME!!!
        :param modulo: Size of the finite field
        """
        self.modulo = modulo

    def get_element(self, int_elem):
        """
        Return a finite field element
        with special operations.
        :param int_elem: Integer value of the element.
        :return: FiniteFieldElement instance
        """
        if not isinstance(int_elem, int):
            raise TypeError(f'Element {int_elem} must be int!')

        return FiniteFieldElement(int_elem, self.modulo)

    def print(self):
        print(f'Finite field of size {self.modulo:,}')

    def get_primitive_elements(self):
        """
        Algorithm to find primitive elements
        of the finite fields.
        Source: https://www.researchgate.net/publication/
        320402737_An_Algorithm_to_Find_Square_Root_of_Quadratic_Residues_over_Finite_Fields_using_Primitive_Elements
        :return: List of primitive elements
        """
        p_elems = []
        for a in range(2, self.modulo):
            i = 0
            q0 = FiniteFieldElement(1, self.modulo)
            while q0.value != 1 or i == 0:
                q0 = q0 * a
                i += 1
            if i == self.modulo - 1:
                p_elems.append(FiniteFieldElement(a, self.modulo))
        return p_elems

    def get_primitive_element(self):
        """
        For square root only one primitive
         element is necessary.
        :return: First primitive element found
        """
        for a in range(2, self.modulo):
            i = 0
            q0 = FiniteFieldElement(1, self.modulo)
            while q0.value != 1 or i == 0:
                q0 = q0 * a
                i += 1
            if i == self.modulo - 1:
                return FiniteFieldElement(a, self.modulo)
