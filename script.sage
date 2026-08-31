from Crypto.Util.number import getPrime
import random
import time

# https://github.com/jvdsn/crypto-attacks/blob/master/shared/polynomial.py
import logging

from sage.all import ZZ
from sage.all import Zmod

# from shared.crt import fast_crt


def _polynomial_hgcd(ring, a0, a1):
    assert a1.degree() < a0.degree()

    if a1.degree() <= a0.degree() / 2:
        return 1, 0, 0, 1

    m = a0.degree() // 2
    b0 = ring(a0.list()[m:])
    b1 = ring(a1.list()[m:])
    R00, R01, R10, R11 = _polynomial_hgcd(ring, b0, b1)
    d = R00 * a0 + R01 * a1
    e = R10 * a0 + R11 * a1
    if e.degree() < m:
        return R00, R01, R10, R11

    q, f = d.quo_rem(e)
    g0 = ring(e.list()[m // 2:])
    g1 = ring(f.list()[m // 2:])
    S00, S01, S10, S11 = _polynomial_hgcd(ring, g0, g1)
    return S01 * R00 + (S00 - q * S01) * R10, S01 * R01 + (S00 - q * S01) * R11, S11 * R00 + (S10 - q * S11) * R10, S11 * R01 + (S10 - q * S11) * R11


def fast_polynomial_gcd(a0, a1):
    """
    Uses a divide-and-conquer algorithm (HGCD) to compute the polynomial gcd.
    More information: Aho A. et al., "The Design and Analysis of Computer Algorithms" (Section 8.9)
    :param a0: the first polynomial
    :param a1: the second polynomial
    :return: the polynomial gcd
    """
    # TODO: implement extended variant of half GCD?
    assert a0.parent() == a1.parent()

    if a0.degree() == a1.degree():
        if a1 == 0:
            return a0
        a0, a1 = a1, a0 % a1
    elif a0.degree() < a1.degree():
        a0, a1 = a1, a0

    assert a0.degree() > a1.degree()
    if a1 == 0:
        return a0.monic()
    ring = a0.parent()

    # Optimize recursive tail call.
    while True:
        logging.debug(f"deg(a0) = {a0.degree()}, deg(a1) = {a1.degree()}")
        _, r = a0.quo_rem(a1)
        if r == 0:
            return a1.monic()

        R00, R01, R10, R11 = _polynomial_hgcd(ring, a0, a1)
        b0 = R00 * a0 + R01 * a1
        b1 = R10 * a0 + R11 * a1
        if b1 == 0:
            return b0.monic()

        _, r = b0.quo_rem(b1)
        if r == 0:
            return b1.monic()

        a0 = b1
        a1 = r
        

def lcg(a, b, x, N):
    return (a*x + b) % N

def gen_params(sz, n):

    p, q = getPrime(sz), getPrime(sz)
    N = p * q

    r = random.randint(1, N)
    Rs = [r]

    a, b = random.randint(1, N), random.randint(1, N)
    for _ in range(n - 1):
        Rs.append(lcg(a, b, Rs[-1], N))

    Cs = [random.randint(1, N) for _ in range(n)]

    s = random.randint(1, N)
    Zs = [ri*pow(s, ci, N) for ri, ci in zip(Rs, Cs)]

    return a, b, s, N, Zs, Cs, Rs

    

def get_poly(N, a, b, combs, Zs):
    Zn.<x> = PolynomialRing(Zmod(N))
    
    lhs = 1
    rhs = 1

    for i in range(len(combs)):
        lhs *= pow(Zs[i], combs[i], N)

    cur = x
    for i in range(len(combs)):
        exp = combs[i]
        if exp < 0:
            lhs *= cur^(-exp)
        else:
            rhs *= cur^exp

        cur = ((cur * a) + b) 

    f = rhs - lhs
    return f
    
    
bit_length = [256, 512, 1024]
sample_size = [200, 150, 100, 50]

for sz in bit_length:
    for n in sample_size:
        
        print('N Bit_length:', 2*sz, 'Sample_size:', n)
        
        for trail_no in range(3):
        
            a, b, s, N, Zs, Cs, Rs = gen_params(sz, n)

            t1 = time.time()

            K = 2^(sz)
            M = matrix(ZZ, [[K * c] for c in Cs]).augment(identity_matrix(len(Cs)))
            L = M.LLL()
            t2 = time.time()
            print('Lattice Reduction time', t2 - t1)


            b_szs = [i for i in L[0][1:]]

            sz1, sz2 = 0, 0

            for bi in b_szs:
                if bi < 0:
                    sz1 += (-bi)
                else:
                    sz2 += bi

            deg = max(sz1, sz2)

            print('poly degree around', deg)

            if deg > 10^7:
                print('TOO BIG')
                print('-'*50)
                continue

            f1 = get_poly(N, a, b, L[0][1:] , Zs)
            t3 = time.time()
            print('Poly find time', t3 - t2)

            print(f1(Rs[0]))
            print(f1.degree())

            f2 = get_poly(N, a, b, L[1][1:], Zs)
            t4 = time.time()
            print('Poly find time', t4 - t3)

            Zn.<x> = PolynomialRing(Zmod(N))

            g = fast_polynomial_gcd(f1, f2)
            t5 = time.time()

            print('GCD time', t5 - t4)


            root = -g[0]/g[1]
            print(root == Rs[0])
            
            r0 = root
            r1 = (a*r0 + b) % N
            
            Zn = Zmod(N)
            v0 = Zs[0] / Zn(r0)
            v1 = Zs[1] / Zn(r1)
            
            _, x, y = xgcd(Cs[0], Cs[1])
            s_ = (Zn(v0)^x) * (Zn(v1)^y)
            print(s == s_)

            print('-'*50)
    
    
