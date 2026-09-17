import sys
from ring import star
from lp import lp, check_css, css_params
from seeds import SEEDS, seed, rep_check, twisted_toric, cluster_cyclic

# reported [[n, k, d]] for each seed, and where it is printed. d is copied as printed (for the
# Cain et al. codes, an upper bound from a numerical search) and is not checked here.
REPORTED = {
    "toy_2x4_l5": ("[[100, 26, 4]]",      "Zheng et al. Sec. III, Eq. (32)"),
    "lp3x5_l33":  ("[[1122, 148, ≤20]]",  "Cain et al. App. A (A3); Zheng et al. Sec. II"),
    "lp3x7_l45":  ("[[2610, 744, ≤16]]",  "Cain et al. App. A (A5)"),
    "lp3x7_l75":  ("[[4350, 1224, ≤20]]", "Cain et al. App. A (A7); Zheng et al. Sec. II"),
    "lp3x7_l91":  ("[[5278, 1480, ≤24]]", "Cain et al. App. A (A9)"),
}

def markdown(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines)

def seed_row(name):
    """One row of the parameter table for LP_l(A, A*)."""
    A = seed(name)
    mA, nA, l = A.shape
    Hx, Hz = lp(A, star(A))
    check_css(Hx, Hz)
    n, k = css_params(Hx, Hz)
    rA = nA - mA  # Zheng et al.'s r_A (Cain et al. use r_A for the number of rows, m_A here)
    kc = rA * rA * l  # logical qubits guaranteed by the construction
    reported, source = REPORTED[name]
    return [name, l, f"{mA}×{nA}", n, k, reported, kc, k - kc,
            Hx.sum(axis=1).max(), Hz.sum(axis=1).max(),
            f"{Hx.sum(axis=0).max()} / {Hz.sum(axis=0).max()}", source]

def sanity_rows():
    """Codes whose parameters are known independently, built before trusting any seed."""
    D3, D7 = rep_check(3, 3), rep_check(7, 33)
    A_cc, B_cc = cluster_cyclic()
    codes = [
        ("twisted toric, lp(A, A)", 3, twisted_toric(3), twisted_toric(3), "[[24, 2]]", "k = 2 for every odd l, Zheng et al. Sec. III"),
        ("twisted toric, lp(A, A)", 5, twisted_toric(5), twisted_toric(5), "[[40, 2]]", "k = 2 for every odd l, Zheng et al. Sec. III"),
        ("stacked surface, lp(D, D*), d=3", 3, D3, star(D3), "[[39, 3]]", "3 copies of the [[13, 1, 3]] surface code"),
        ("stacked surface, lp(D, D*), d=7", 33, D7, star(D7), "[[2805, 33, 7]]", "Zheng et al. Table V"),
        ("cluster cyclic, lp(A, B)", 5, A_cc, B_cc, "[[40, 8, 5]]", "Zheng et al. Eq. (20)"),
    ]
    rows = []
    for label, l, A, B, expected, source in codes:
        Hx, Hz = lp(A, B)
        check_css(Hx, Hz)
        n, k = css_params(Hx, Hz)
        rows.append([label, l, n, k, expected, source])
    return rows

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")  # the table uses ≤, ×, ² and −
    print(markdown(["code", "l", "seed A", "n", "k", "reported", "r_A² l", "k − r_A² l",
                    "max X-check weight", "max Z-check weight", "max qubit degree X / Z", "reported in"],
                   [seed_row(name) for name in SEEDS]))
    print()
    print(markdown(["sanity code", "l", "n", "k", "expected", "source"], sanity_rows()))
