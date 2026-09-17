This repository reproduces the lifted product (LP) codes from the paper "Logical computation with canonical lifted product codes" by Han Zheng, Guo Zheng, Liang Jiang and Qian Xu ([arXiv:2607.28605](https://arxiv.org/abs/2607.28605)), built from their definitions and checked against the stated parameters. That paper does not print the seed matrices of its codes; they are in its reference [41], "Shor’s algorithm is possible with as few as 10,000 reconfigurable atomic qubits" by Madelyn Cain, Qian Xu and others ([arXiv:2603.28627](https://arxiv.org/abs/2603.28627)), Appendix A.

## Results

| code | l | seed A | n | k | reported | r_A² l | k − r_A² l | max X-check weight | max Z-check weight | max qubit degree X / Z | reported in |
|---|---|---|---|---|---|---|---|---|---|---|---|
| toy_2x4_l5 | 5 | 2×4 | 100 | 26 | [[100, 26, 4]] | 20 | 6 | 6 | 6 | 4 / 4 | Zheng et al. Sec. III, Eq. (32) |
| lp3x5_l33 | 33 | 3×5 | 1122 | 148 | [[1122, 148, ≤20]] | 132 | 16 | 8 | 8 | 5 / 5 | Cain et al. App. A (A3); Zheng et al. Sec. II |
| lp3x7_l45 | 45 | 3×7 | 2610 | 744 | [[2610, 744, ≤16]] | 720 | 24 | 10 | 10 | 7 / 7 | Cain et al. App. A (A5) |
| lp3x7_l75 | 75 | 3×7 | 4350 | 1224 | [[4350, 1224, ≤20]] | 1200 | 24 | 10 | 10 | 7 / 7 | Cain et al. App. A (A7); Zheng et al. Sec. II |
| lp3x7_l91 | 91 | 3×7 | 5278 | 1480 | [[5278, 1480, ≤24]] | 1456 | 24 | 10 | 10 | 7 / 7 | Cain et al. App. A (A9) |

For all five codes, my $n$ and $k$ equal the reported values, and the CSS condition $H_X H_Z^T = 0 \pmod 2$ holds. The check weights (8 for the 3×5 seed, 10 for the 3×7 seeds) are the stabilizer weights $3+5$ and $3+7$ given by Cain et al., and the qubit degrees match Zheng et al. Tables II and V. None of these depend on $l$.

**Distances are not checked.** The $d$ in the "reported" column is copied from the papers. For the Cain et al. codes it is an upper bound from a numerical search.

In every code, $k$ is larger than the $r_A^2 l$ logical qubits the construction guarantees: by 6, 16, 24, 24 and 24. The excess is the same 24 for three different 3×7 seeds at three different lift sizes. Cain et al. only remark that $k$ exceeds the bound. For the toy code, Zheng et al. attribute the 6 extra logical qubits to the $1+x$ factor of $x^5+1$. I have not yet worked out where 16 and 24 come from.