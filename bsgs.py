from elliptic_curve import ECPoint, EllipticCurve
from finite_field import FiniteFieldElement, FiniteField
from math import ceil, sqrt
import time


def generate_baby_steps(p, m):
    """
    Generate baby step list,
    of multiples of P: a*P.
    :param p: ECPoint P
    :param m: number of babysteps
    :return:
    """
    print(f'Generating {m} baby steps...')
    return [a*p for a in range(m+1)]


def giant_steps(p, q, baby_steps, m):
    """
    Find collision of BS and calculated GS.
    :param p: ECPoint P
    :param q: ECPoint Q
    :param baby_steps: list of pregenerated baby steps
    :param m: number of babysteps
    :return: Result of the algorithm - log_P Q
    """
    j = 0
    i = 0
    while True:
        x = q - j*p
        try:
            i = baby_steps.index(x)
        except ValueError:
            j += 1
            continue
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
    r = p.order()
    m = ceil(sqrt(r))

    begin = time.time()
    baby_steps = generate_baby_steps(p, m)
    end = time.time()

    print(f'Babysteps generated in {(end-begin):.3f} seconds.\n')

    p2 = m*p

    res_list = []

    for i in q_list:
        begin = time.time()
        res_list.append(giant_steps(p2, i, baby_steps, m))
        end = time.time()
        print(f'Logarithm of point {i} found in {(end-begin):.3f} seconds.\n')

    return res_list
