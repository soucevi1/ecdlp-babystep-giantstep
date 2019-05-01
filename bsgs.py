from elliptic_curve import ECPointAtInfinity, binary_search
from math import ceil, sqrt
import time


class BabyStepPoint:

    def __init__(self, point, index):
        self.point = point
        self.index = index

    def __lt__(self, other):
        return self.point < other.point


def generate_baby_steps(p, m):
    """
    Generate baby step list,
    of multiples of P: a*P.
    :param p: ECPoint P
    :param m: number of babysteps
    :return:
    """
    print(f'Generating {m:_} baby steps...')
    # return [a*p for a in range(m-1)]
    b = ECPointAtInfinity(p.curve)
    baby_steps = [BabyStepPoint(b, 0)]
    for i in range(1, m-1):
        b = p+b
        baby_steps.append(BabyStepPoint(b, i))
    print('Sorting baby steps for more efficient search...')
    baby_steps.sort()
    return baby_steps


def giant_steps(p, q, baby_steps, m):
    """
    Find collision of BS and calculated GS.
    :param p: ECPoint P
    :param q: ECPoint Q
    :param baby_steps: list of pre-generated baby steps
    :param m: number of babysteps
    :return: Result of the algorithm - log_P Q
    """
    j = 0
    new_p = ECPointAtInfinity(p.curve)
    while True:
        # x = q - j*p
        x = q - new_p
        i = binary_search(baby_steps, x)
        new_p = new_p + p
        if i == -1:
            # No collision found
            j += 1
            continue
        print(f'Collision! Index in baby steps: {i}, index in giant steps: {j}')
        result = i + j*m
        return result


def find_logarithm(q_list, p):
    """
    Find n such that n*P = Q for each Q
    in the q_list
    :param q_list: list of ECPoints Q
    :param p: ECPoint P
    :return: list of logarithms for all Qs
    """
    # r = p.order()
    r = p.order_approx()
    m = ceil(sqrt(r))

    begin = time.time()
    baby_steps = generate_baby_steps(p, m)
    end = time.time()

    print(f'Babysteps generated and sorted in {(end-begin):.3f} seconds.\n')

    p2 = m*p

    res_list = []

    for i in q_list:
        begin = time.time()
        res_list.append(giant_steps(p2, i, baby_steps, m))
        end = time.time()
        print(f'Logarithm of point {i} found in {(end-begin):.3f} seconds.\n')

    return res_list
