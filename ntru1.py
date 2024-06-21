from fractions import Fraction as frac
from fractions import Fraction as frac
from operator import add
from operator import neg
from operator import mod

# Helper Functions
def modPoly(c, k):
    if k == 0:
        raise ValueError("Integer k must be non-zero")
    return [fracMod(x, k) for x in c]

def subPoly(c1, c2):
    [c1, c2] = resize(c1, c2)
    c2 = list(map(neg, c2))
    out = list(map(add, c1, c2))
    return trim(out)

def multPoly(c1, c2):
    order = (len(c1) - 1 + len(c2) - 1)
    out = [0] * (order + 1)
    for i in range(len(c1)):
        for j in range(len(c2)):
            out[j + i] += c1[i] * c2[j]
    return trim(out)

def resize(c1, c2):
    if len(c1) > len(c2):
        c2 = c2 + [0] * (len(c1) - len(c2))
    if len(c1) < len(c2):
        c1 = c1 + [0] * (len(c2) - len(c1))
    return [c1, c2]

def trim(seq):
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] != 0:
            break
    return seq[:i + 1]

def extEuclidPoly(a, b):
    switch = False
    a = trim(a)
    b = trim(b)
    if len(a) >= len(b):
        a1, b1 = a, b
    else:
        a1, b1 = b, a
        switch = True
    Q, R = [], []
    while b1 != [0]:
        [q, r] = divPoly(a1, b1)
        Q.append(q)
        R.append(r)
        a1 = b1
        b1 = r
    S = [0] * (len(Q) + 2)
    T = [0] * (len(Q) + 2)
    S[0], S[1], T[0], T[1] = [1], [0], [0], [1]
    for x in range(2, len(S)):
        S[x] = subPoly(S[x-2], multPoly(Q[x-2], S[x-1]))
        T[x] = subPoly(T[x-2], multPoly(Q[x-2], T[x-1]))
    gcdVal = R[len(R)-2]
    s_out = S[len(S)-2]
    t_out = T[len(T)-2]
    scaleFactor = gcdVal[len(gcdVal)-1]
    gcdVal = [x/scaleFactor for x in gcdVal]
    s_out = [x/scaleFactor for x in s_out]
    t_out = [x/scaleFactor for x in t_out]
    if switch:
        return [gcdVal, t_out, s_out]
    else:
        return [gcdVal, s_out, t_out]

def divPoly(N, D):
    N, D = list(map(frac, trim(N))), list(map(frac, trim(D)))
    degN, degD = len(N) - 1, len(D) - 1
    if degN >= degD:
        q = [0] * (degN - degD + 1)
        while degN >= degD and N != [0]:
            d = list(D)
            [d.insert(0, frac(0, 1)) for i in range(degN - degD)]
            q[degN - degD] = N[degN] / d[len(d) - 1]
            d = [x * q[degN - degD] for x in d]
            N = subPoly(N, d)
            degN = len(N) - 1
        r = N
    else:
        q = [0]
        r = N
    return [trim(q), trim(r)]

def addPoly(c1, c2):
    [c1, c2] = resize(c1, c2)
    out = list(map(add, c1, c2))
    return trim(out)

def cenPoly(c, q):
    u = float(q) / float(2)
    l = -u
    c = modPoly(c, q)
    c = [x - q if x > u else x for x in c]
    c = [x + q if x <= l else x for x in c]
    return c

def reModulo(num, div, modby):
    [_, remain] = divPoly(num, div)
    return modPoly(remain, modby)

def egcd(a, b):
    x, y, u, v = 0, 1, 1, 0
    while a != 0:
        q, r = b // a, b % a
        m, n = x - u * q, y - v * q
        b, a, x, y, u, v = a, r, u, v, m, n
    gcdVal = b
    return gcdVal, x, y

def modinv(a, m):
    gcdVal, x, y = egcd(a, m)
    if gcdVal != 1:
        return None  # modular inverse does not exist
    else:
        return x % m

def fracMod(f, m):
    [tmp, _, _] = egcd(f.denominator, m)
    if tmp != 1:
        raise ValueError("GCD of denominator and m is not 1")
    else:
        out = modinv(f.denominator, m) * f.numerator % m
        return out

# NTRU Parameters
N = 11
p = 3
q = 32

# Simpler polynomials
f = [1, -1, 1] #f = 1-x+x^2
g = [1, 1]     #g = 1+x

D = [0] * (N + 1)
D[0] = -1
D[N] = 1

print("Key Generation")
print("Values used:")
print(" N =", N)
print(" p =", p)
print(" q =", q)
print("\n")
print("two polynomials (g and f):")
print("f(x)= ", f)
print("g(x)= ", g)

print("\nNow we determine Fp and Fq")
[gcd_f, s_f, t_f] = extEuclidPoly(f, D)

f_p = modPoly(s_f, p)
f_q = modPoly(s_f, q)
print("Fp:", f_p)
print("Fq:", f_q)

x = multPoly(f_q, g)
h = reModulo(x, D, q)

print("\nh is determined")
print("fq x g: ", x)
print("H (Public Key): ", h)

print("\nEncryption")
msg = [1, 0, 1,1,0,1,1]
randPol = [0, 1,1,-1,1,1,1]

print("Message:\t\t", msg)
print("Random:\t\t\t", randPol)
e_tilda = addPoly(multPoly(multPoly([p], randPol), h), msg)
e = reModulo(e_tilda, D, q)

print("Encrypted message:\t", e)

print("\nDecryption")

tmp = reModulo(multPoly(f, e), D, q)
centered = cenPoly(tmp, q)
m1 = multPoly(f_p, centered)
tmp = reModulo(m1, D, p)

print("Decrypted message:\t", trim(tmp))