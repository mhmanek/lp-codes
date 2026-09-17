import numpy as np

def row_reduce(A):
    """Reduced row echelon form of A over GF(2).
    Returns (R, pivots) where pivots is the list of pivot column indices."""
    A = A.copy().astype(np.uint8) % 2 # first convert into a binary matrix by taking mod 2
    r = 0 # initialise variable for pivot-row counter r
    pivots = []
    for c in range(np.shape(A)[1]): # walk left to right over columns
        for j in range(r, np.shape(A)[0]): # for each column: find a row at index >= r with a 1 there; if none, skip the column
            if A[j, c] != 0:
                A[[r, j]] = A[[j, r]] # swap that row up to position r
                for i in range(0, np.shape(A)[0]):
                    if i != r and A[i, c] != 0:
                        A[i] ^= A[r] # XOR row r into every other row that has a 1 in this column
                pivots.append(c) # append the column index to pivots
                r += 1
                break
    return(A, pivots) 

def rank(A):
    """Rank of A over GF(2). Each row is packed into one Python integer, so XOR-ing two rows is a
    single operation: a few seconds for a 1911 x 5278 check matrix instead of minutes with row_reduce."""
    A = np.asarray(A, dtype=np.uint8) % 2
    basis = {} # leading bit position -> the stored row whose highest set bit is at that position
    for row in A:
        v = int.from_bytes(np.packbits(row).tobytes(), "big") # the row as one big integer
        while v:
            lead = v.bit_length() - 1 # position of the highest set bit
            if lead not in basis:
                basis[lead] = v # new leading bit: this row is independent of the rows stored so far
                break
            v ^= basis[lead] # clear that bit and keep reducing
    return len(basis) # the rank is how many independent rows were stored

def nullspace(A):
    """Basis for {x : A x = 0 mod 2}, returned as rows of a matrix."""
    R, pivots = row_reduce(A)
    n = A.shape[1]
    pivot_set = set(pivots)
    free_variables = [c for c in range(n) if c not in pivot_set]
    basis = []
    for j in free_variables:
        v = np.zeros(n, dtype=np.uint8)
        v[j] = 1
        for k in range(len(pivots)):
            v[pivots[k]] = R[k, j]
        basis.append(v)
    return np.array(basis, dtype=np.uint8).reshape(len(basis), n)

def in_rowspace(A, v):
    """Is the vector v a linear combination of the rows of A?"""
    return rank(A) == rank(np.vstack([A, v])) # compare rank(A) with rank of (A with v stacked underneath it)

def solve(M, b):
    """One solution g of M g = b over GF(2), or None if there is none."""
    M = np.asarray(M, dtype=np.uint8) % 2
    b = np.asarray(b, dtype=np.uint8).reshape(-1, 1) % 2
    n = M.shape[1]
    R, pivots = row_reduce(np.hstack([M, b]))
    if n in pivots:
        return None
    g = np.zeros(n, dtype=np.uint8)
    for k, c in enumerate(pivots):
        g[c] = R[k, n]
    return g