from ring import kron, ident, star, lift
from f2linalg import rank
import numpy as np

def lp(A, B):
    """Lifted product LP_l(A, B), paper's Eq. (1). Returns binary (Hx, Hz)."""
    mA, nA, l = A.shape
    mB, nB, _ = B.shape
    HX = np.concatenate([kron(A, ident(mB, l)), kron(ident(mA, l), B)], axis=1)
    HZ = np.concatenate([kron(ident(nA, l), star(B)), kron(star(A), ident(nB, l))], axis=1)
    return lift(HX), lift(HZ)

def qubit_index(block, i, j, m, A, B):
    """Flat column index of qubit q_L(i,j,m) (block=0) or q_R(i,j,m) (block=1)."""
    mA, nA, l = A.shape
    mB, nB, _ = B.shape
    if block == 0:
        # left block: fibre (i, j) with i in [nA], j in [mB]
        return (i * mB + j) * l + m
    # right block starts after the nA*mB*l left-block qubits
    # fibre (i, j) with i in [mA], j in [nB]
    return nA * mB * l + (i * nB + j) * l + m

def css_params(Hx, Hz):
    n = Hx.shape[1]
    assert n == Hz.shape[1]
    k = n - rank(Hx) - rank(Hz)
    return n, k

def check_css(Hx, Hz):
    assert Hx.shape[1] == Hz.shape[1]
    # float64 so numpy hands the matmul to BLAS: int64 matmul has no fast path (16 s vs 0.3 s on the
    # 1911 x 5278 matrices). Entries are overlap counts of at most n, exact in float64 (and no uint8 overflow).
    product = (Hx.astype(np.float64) @ Hz.T.astype(np.float64)) % 2
    assert not product.any(), "CSS condition fails: Hx Hz^T is not 0 mod 2"