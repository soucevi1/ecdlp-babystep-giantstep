from math import gcd


class BabyStepPoint:
    """
    Helper class to store generated
    babysteps in a sorted list and
    keep their indexes after sorting.
    """

    def __init__(self, point, index):
        self.point = point
        self.index = index

    def __lt__(self, other):
        return self.point < other.point


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
    return abs(a * b) // gcd(a, b)


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
        mid = (lo + hi) // 2
        midval = element_list[mid].point
        if midval < element:
            lo = mid + 1
        elif midval > element:
            hi = mid
        else:
            return element_list[mid].index
    return -1
