import bsgs
from finite_field import FiniteField, FiniteFieldElement
from elliptic_curve import EllipticCurve, ECPoint


def main():
    p = 7919

    print('Creating finite field: ')
    ff = FiniteField(p)
    ff.print()
    print('')

    a = ff.get_element(0)
    b = ff.get_element(1)

    print('Creating elliptic curve: ')
    e = EllipticCurve(a, b, ff)
    print(e)
    print('')

    P = ECPoint(6692, 191, e)
    print(f'Point P: {P}')

    m1 = 496
    m2 = 1759

    Q1 = m1*P
    Q2 = m2*P

    qlist = [Q1, Q2]
    for i in range(len(qlist)):
        print(f'Point Q{i+1}: {qlist[i]}')
    print('')

    results = bsgs.find_logarithm(qlist, P)

    for i in range(len(results)):
        print('================================================')
        print(f'lop_P Q{i+1} = {results[i]}')
        print(f'Result is correct: {results[i] * P == qlist[i]}')

    print('================================================')


if __name__ == '__main__':
    main()
