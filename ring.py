import numpy as np

def poly(exponents, l):
    """Ring element with a 1 at each listed exponent, reduced mod l.
    poly([0, 3], 5) is 1 + x^3, stored as [1 0 0 1 0]."""
    p = np.zeros(l, dtype=np.uint8)
    for e in exponents:
        p[e % l] ^= 1
    return p

def pstr(a):
    """Readable form of a ring element: [1 0 0 1 0] -> '1 + x^3'.
    The zero element gives '0'."""
    terms = []
    for e in range(len(a)):
        if a[e] == 1:
            if e == 0:
                terms.append("1")
            elif e == 1:
                terms.append("x")
            else:
                terms.append("x^" + str(e))
    if len(terms) == 0:
        return "0"
    return " + ".join(terms)

def pmul(a, b):
    """Product in R: cyclic convolution mod 2."""
    l = len(a)
    result = np.zeros(l, dtype=np.uint8)
    for s in range(l):
        if a[s] == 1:
            result ^= np.roll(b, s)
    return result

def conj(a):
    """a(x^-1): the coefficient at index e moves to index (-e) % l."""
    l = len(a)
    out = np.zeros(l, dtype=np.uint8)
    for e in range(l):
        out[(-e) % l] = a[e]
    return out

def circ(a):
    """The l x l binary matrix B(a) of multiplication by a.
    Entry [r, t] is a[(r - t) % l]."""
    l = len(a)
    C = np.zeros((l, l), dtype=np.uint8)
    for r in range(l):
        for t in range(l):
            C[r, t] = a[(r - t) % l]
    return C

def ring_matrix(exponent_table, l):
    """exponent_table[i][j] is a list of exponents ([] for a zero entry).
    ring_matrix([[[0], [1]], [[0, 1], []]], 3) is [[1, x], [1+x, 0]],
    of shape (2, 2, 3)."""
    rows = len(exponent_table)
    cols = len(exponent_table[0])
    M = np.zeros((rows, cols, l), dtype=np.uint8)
    for i in range(rows):
        for j in range(cols):
            M[i, j] = poly(exponent_table[i][j], l)
    return M

def ident(n, l):
    """n x n matrix over R, the constant polynomial 1 on the diagonal.
    Shape (n, n, l)."""
    M = np.zeros((n, n, l), dtype=np.uint8)
    for i in range(n):
        M[i, i, 0] = 1
    return M

def star(M):
    """M*: swap the (rows, cols) axes and conj every entry.
    Shape (m, n, l) -> (n, m, l)."""
    m, n, l = M.shape
    out = np.zeros((n, m, l), dtype=np.uint8)
    for i in range(m):
        for j in range(n):
            out[j, i] = conj(M[i, j])      # not M.T: that would reverse all 3 axes
    return out

def kron(A, B):
    """Kronecker product over R:
    entry [i*rB + k, j*cB + m] is pmul(A[i, j], B[k, m])."""
    rA, cA, l = A.shape
    rB, cB, _ = B.shape
    out = np.zeros((rA * rB, cA * cB, l), dtype=np.uint8)
    for i in range(rA):
        for j in range(cA):
            if not A[i, j].any():
                continue
            for k in range(rB):
                for m in range(cB):
                    if not B[k, m].any():
                        continue
                    out[i * rB + k, j * cB + m] = pmul(A[i, j], B[k, m])
    return out

def lift(M):
    """Binarization B(M): replace every entry by its circulant.
    Shape (rows*l, cols*l)."""
    rows, cols, l = M.shape
    out = np.zeros((rows * l, cols * l), dtype=np.uint8)
    for i in range(rows):
        for j in range(cols):
            out[i * l:(i + 1) * l, j * l:(j + 1) * l] = circ(M[i, j])
    return out