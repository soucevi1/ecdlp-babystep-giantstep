
# The main module of the ECDLP-BSGS program.
# Author: Vit Soucek

import bsgs
from elliptic_curve import EllipticCurve, ECPoint
from finite_field import FiniteField


def test_initialize():
    """
    Initialize the logarithm parameters
    with small values for testing purposes.
    :return: Value of P, Q1 and Q2, where we search log_P Q*
    """
    finite_field_prime = 7919

    print('Creating finite field: ')
    ff = FiniteField(finite_field_prime)
    ff.print()
    print('')

    print('Creating elliptic curve: ')
    a = ff.get_element(0)
    b = ff.get_element(1)
    e = EllipticCurve(a, b, ff)
    print(e)
    print('')

    point_p = ECPoint(6692, 191, e)
    print(f'Point P: {point_p}')

    m1 = 496
    m2 = 1759

    point_q1 = m1 * point_p
    point_q2 = m2 * point_p

    return point_p, point_q1, point_q2


def initialize():
    """
    Initialize the logarithm parameters
    with actual values given in the task.
    :return: Value of P, Q1 and Q2, where we search log_P Q*
    """
    finite_field_prime = 2 ** 42 + 1597

    print('Creating finite field: ')
    ff = FiniteField(finite_field_prime)
    ff.print()
    print('')

    print('Creating elliptic curve: ')
    a = ff.get_element(0)
    b = ff.get_element(1)
    e = EllipticCurve(a, b, ff)
    print(e)
    print('')

    point_p = ECPoint(3, 678235393584, e)
    print(f'Point P: {point_p}')

    point_q1 = ECPoint(3480506547504, 2007642106263, e)
    point_q2 = ECPoint(3563933383778, 1990746467823, e)

    return point_p, point_q1, point_q2


def main():
    print(f'Searching x such as: x = log_P Q')
    print(f'--------------------------------')
    point_p, point_q1, point_q2 = initialize()
    # point_p, point_q1, point_q2 = test_initialize()

    q_list = [point_q1, point_q2]

    for i in range(len(q_list)):
        print(f'Point Q{i + 1}: {q_list[i]}')
    print('')

    results = bsgs.find_logarithm(q_list, point_p)

    for i in range(len(results)):
        print('================================================')
        print(f'lop_P Q{i + 1} = {results[i]}')
        print(f'Result is correct: {results[i] * point_p == q_list[i]}')

    print('================================================')


if __name__ == '__main__':
    main()
