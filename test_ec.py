from finite_field import FiniteField, FiniteFieldElement
from elliptic_curve import EllipticCurve, ECPoint, ECPointAtInfinity

f = FiniteField(7919)

print('FF: ')
f.print()

a = f.get_element(0)
b = f.get_element(1)

e = EllipticCurve(a, b, f)

p = ECPoint(f.get_element(6692), f.get_element(191), e)
q = ECPoint(f.get_element(4445), f.get_element(5180), e)

print(f'p+q: {p+q}')
print(f'5*p: {5*p}')

print(f'order of p: {p.order()}')
print(f'order of q: {q.order()}')

p = e.random_point()
print(f'random point: {p}')
print(f'is it a point? {e.is_point(p.x, p.y)}')
