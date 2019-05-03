
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


def binary_search(element_list, element):
    """
    Binary search algorithm in list of objects.
    Source: https://stackoverflow.com/a/4161779/6136143
    :param element_list: List of elements (BabyStepPoints)
    :param element: Element to look for (ECPoint)
    :return: Index of element, or -1 in case of unsuccessful search.
    """
    hi = len(element_list)
    lo = 0
    while lo < hi:
        mid = (lo + hi) // 2
        mid_val = element_list[mid].point
        if mid_val < element:
            lo = mid + 1
        elif mid_val > element:
            hi = mid
        else:
            return element_list[mid].index
    return -1
