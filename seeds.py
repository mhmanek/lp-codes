from ring import ring_matrix, poly

# each seed: (l, exponent table). x^0 = 1. Source equation in the comment.
SEEDS = {
    "toy_2x4_l5":   (5,  [[0,0,0,0], [0,1,2,4]]),                # Zheng et al. Eq. (32)
    "lp3x5_l33":    (33, [[0,0,0,0,0], [0,14,19,11,26], [0,13,2,15,21]]),   # Cain et al. (A3)
    "lp3x7_l45":  (45, [[29,21,31,15,37,25,27], [13,25,19,26,11,18,29], [31,2,27,32,41,41,18]]), # Cain et al. (A5)
    "lp3x7_l75":  (75, [[0,71,73,68,33,50,47], [38,39,60,26,18,1,23], [73,6,5,42,20,22,73]]), # Cain et al. (A7)
    "lp3x7_l91":  (91, [[57,75,42,80,7,67,27], [57,73,34,12,27,50,87], [21,53,70,18,1,3,18]]) # Cain et al. (A9)
}

def seed(name):
    """Return the ring matrix A for a named seed (every entry a single monomial)."""
    l, table = SEEDS[name]
    return ring_matrix([[[e] for e in row] for row in table], l)

def rep_check(d, l):
    """(d-1) x d repetition-code check matrix over R, entries the constant 1."""
    # row i has the constant 1 in columns i and i+1, zero elsewhere
    table = [[[0] if j in (i, i + 1) else [] for j in range(d)] for i in range(d - 1)]
    return ring_matrix(table, l)

def twisted_toric(l):
    # [[1, 1], [1, x]]
    return ring_matrix([[[0], [0]], [[0], [1]]], l)

def cluster_cyclic():
    """(A, B) of Zheng et al. Eq. (20), l = 5. B is not A*, so this checks lp(A, B) in general."""
    # A = [[x+x^3, x+x^4], [x+x^4, x+x^4]],  B = [[1+x, x^3+x^4], [1+x^4, 1+x^4]]
    A = ring_matrix([[[1, 3], [1, 4]], [[1, 4], [1, 4]]], 5)
    B = ring_matrix([[[0, 1], [3, 4]], [[0, 4], [0, 4]]], 5)
    return A, B
