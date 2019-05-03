
# Module for babystep-giantstep ECDLP calculation.
# Author: Vit Soucek

from elliptic_curve import ECPointAtInfinity
from math import ceil, sqrt
import time
from helper_tools import BabyStepPoint, binary_search


def generate_baby_steps(p, m):
    """
    Generate baby step list,
    of multiples of P: a*P.
    :param p: ECPoint P
    :param m: number of babysteps
    :return:
    """
    print(f'Generating {m:,} baby steps...')

    b = ECPointAtInfinity(p.curve)
    baby_steps = [BabyStepPoint(b, 0)]

    for i in range(1, m - 1):
        b = p + b
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

        if i != -1:
            print(f'Collision! Index in baby steps: {i}, index in giant steps: {j}')
            result = i + j * m
            return result

        # No collision found
        j += 1


def find_logarithm(q_list, p):
    """
    Find n such that n*P = Q for each Q
    in the q_list
    :param q_list: list of ECPoints Q
    :param p: ECPoint P
    :return: list of logarithms for all Qs
    """
    r = p.order_approx()
    m = ceil(sqrt(r))

    begin = time.time()
    baby_steps = generate_baby_steps(p, m)
    end = time.time()

    print(f'Babysteps generated and sorted in {(end - begin):.3f} seconds.\n')

    p2 = m * p

    res_list = []

    for i in q_list:
        begin = time.time()
        res_list.append(giant_steps(p2, i, baby_steps, m))
        end = time.time()
        print(f'Logarithm of point {i} found in {(end - begin):.3f} seconds.\n')

    return res_list
