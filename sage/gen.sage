from random import randrange

load("field.sage")

# Define curve coefficients, notice we are using a different b than the original paper.
q = ZZ(2^127)
one = F2m(1)
a = F2x(s);

for i in range(0, 127):
	b = F2m(z^i)
	E = EllipticCurve(F2m, [1, 1, 0, 0, b])
	n = E.order();
	t = q + 1 - n
	n2 = (q - 1)^2 + t^2
	r = n2//2
	E2 = EllipticCurve(F2x, [1, s, 0, 0, b])
	if (r in Primes()):
		print(E2)
		break;

if (r in Primes()):
	print(E)
	exit(0)

for i in range(1, 127):
	b = F2m(z^i + 1)
	E = EllipticCurve(F2m, [1, 1, 0, 0, b])
	n = E.cardinality();
	t = q + 1 - n
	n2 = (q - 1)^2 + t^2
	E2 = EllipticCurve(F2x, [1, s, 0, 0, b])
	r = n2//2
	if (r in Primes()):
		print(E2)
		exit(0);


