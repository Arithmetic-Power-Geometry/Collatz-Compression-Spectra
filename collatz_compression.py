"""Core routines for compression-spectrum analysis of Collatz/Syracuse dynamics.

All integer dynamics use Python arbitrary-precision integers.  The main map is the
accelerated odd (Syracuse) map
    S(n) = (3 n + 1) / 2^{nu_2(3 n + 1)}.
The module also provides the exact compression/residual decomposition and sharp
fixed-total-valuation residual bounds used in the accompanying manuscript.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import log
import math
from fractions import Fraction
import numpy as np

LN2 = log(2.0)
LN3 = log(3.0)


def v2(x: int) -> int:
    """2-adic valuation of a positive integer."""
    if x <= 0:
        raise ValueError("v2 expects a positive integer")
    return (x & -x).bit_length() - 1


def odd_core(n: int) -> int:
    if n <= 0:
        raise ValueError("n must be positive")
    return n >> v2(n)


def accelerated_step(n: int, a: int = 3, b: int = 1) -> tuple[int, int]:
    if n <= 0 or n % 2 == 0:
        raise ValueError("accelerated_step expects a positive odd n")
    if a % 2 == 0 or b % 2 == 0:
        raise ValueError("a and b must be odd")
    y = a * n + b
    k = v2(y)
    return y >> k, k


def accelerated_block(n: int, m: int, a: int = 3, b: int = 1):
    """Return values, valuations, cumulative compression and residuals."""
    if m < 1:
        raise ValueError("m must be >= 1")
    x = odd_core(n)
    n0 = x
    vals = [x]
    ks, cs, rs = [], [], []
    K = 0
    lna = log(float(a))
    for r in range(1, m + 1):
        x, k = accelerated_step(x, a, b)
        K += k
        c = K * LN2 - r * lna
        rr = log(x / n0) + c
        vals.append(x); ks.append(k); cs.append(c); rs.append(rr)
    return vals, np.asarray(ks, dtype=np.int16), np.asarray(cs), np.asarray(rs)


def collatz_raw_metrics(n: int, max_steps: int = 10000) -> tuple[int, int, bool]:
    x = int(n); mx = x; steps = 0
    while x != 1 and steps < max_steps:
        x = 3 * x + 1 if (x & 1) else x >> 1
        if x > mx: mx = x
        steps += 1
    return steps, mx, x == 1


def iid_anticompression_probability(m: int) -> float:
    """Exact iid geometric(1/2) baseline P(C_m < 0)."""
    cutoff = math.ceil(m * math.log(3, 2)) - 1
    if cutoff < m:
        return 0.0
    return float(sum(math.comb(K - 1, m - 1) * 2.0 ** (-K)
                     for K in range(m, cutoff + 1)))


def geometric_cramer_rate_mean_k(x: float) -> float:
    """Cramer rate J(x) for the mean of geometric(1/2) variables on {1,2,...}."""
    if x < 1:
        return math.inf
    if abs(x - 1.0) < 1e-15:
        return LN2
    return x * math.log(2.0 * (x - 1.0) / x) - math.log(x - 1.0)


def geometric_cramer_rate_lambda(lam: float) -> float:
    """Cramer rate for Lambda = (mean(k))*log 2 - log 3."""
    x = (lam + LN3) / LN2
    return geometric_cramer_rate_mean_k(x)


def residual_kernel_from_word(word: list[int] | tuple[int, ...]) -> float:
    """S(k)=sum_{i=0}^{m-1} 2^{K_i}/3^i entering R=log(1+S/(3n0))."""
    K = 0
    s = 1.0
    for i, k in enumerate(word[:-1], start=1):
        K += int(k)
        s += (2.0 ** K) / (3.0 ** i)
    return s


def residual_kernel_extrema(m: int, K: int) -> tuple[float, float]:
    """Sharp min/max of residual kernel over positive valuation words of length m, total K.

    For k_j >= 1 and sum k_j = K >= m:
      min: (1,...,1,K-m+1)  (all excess back-loaded)
      max: (K-m+1,1,...,1)  (all excess front-loaded)
    """
    if m < 1 or K < m:
        raise ValueError("require m >= 1 and K >= m")
    smin = 3.0 * (1.0 - (2.0 / 3.0) ** m)
    if m == 1:
        smax = 1.0
    else:
        smax = 1.0 + (2.0 ** (K - m + 1)) * (1.0 - (2.0 / 3.0) ** (m - 1))
    return smin, smax


def residual_bounds_fixed_total(n0: int, m: int, K: int) -> tuple[float, float]:
    """Sharp word-uniform residual bounds for fixed (m,K)."""
    smin, smax = residual_kernel_extrema(m, K)
    return math.log1p(smin / (3.0 * n0)), math.log1p(smax / (3.0 * n0))


def uniform_descent_certificate(n0: int, m: int, K: int) -> bool:
    """Sufficient certificate: every positive valuation word with (m,K) descends by time m."""
    _, rmax = residual_bounds_fixed_total(n0, m, K)
    C = K * LN2 - m * LN3
    return C > rmax




def near_extremal_count_bound(m: int, K: int, epsilon: float) -> dict:
    """Rigorous upper bound for near-maximal residual-kernel words.

    A word is epsilon-near-extremal when
        S(word) >= (1-epsilon) S_max(m,K).
    Writing L=K-m and q for excess valuation delayed beyond the first
    coordinate, Theorem 15 of the manuscript gives
        S_max-S >= 2^(L+1)/3 * (1-2^(-q)).
    If alpha < 1, this bounds q by Q and hence bounds the number of
    near-extremal words by C(Q+m-1,m-1).  When alpha >= 1 the routine
    returns the trivial all-compositions bound.
    """
    if m < 2 or K < m:
        raise ValueError("require m >= 2 and K >= m")
    if not (0.0 <= epsilon < 1.0):
        raise ValueError("epsilon must lie in [0,1)")
    L = K - m
    total = math.comb(K - 1, m - 1)
    if L == 0:
        return {
            'm': m, 'K': K, 'epsilon': epsilon, 'alpha': 0.0, 'Q': 0,
            'count_bound': 1, 'total_words': total, 'fraction_bound': 1.0,
            'nontrivial': True
        }
    _, smax = residual_kernel_extrema(m, K)
    alpha = epsilon * 3.0 * smax / (2.0 ** (L + 1))
    if alpha >= 1.0:
        Q = L
        count_bound = total
        nontrivial = False
    else:
        Q = min(L, math.floor(-math.log2(1.0 - alpha)))
        count_bound = math.comb(Q + m - 1, m - 1)
        nontrivial = True
    return {
        'm': m, 'K': K, 'epsilon': epsilon, 'alpha': alpha, 'Q': Q,
        'count_bound': count_bound, 'total_words': total,
        'fraction_bound': count_bound / total, 'nontrivial': nontrivial
    }

def valuation_word_affine_constant(word: list[int] | tuple[int, ...]) -> tuple[int, int]:
    """Return the exact affine constant A and total valuation K for a word.

    After m accelerated odd steps,
        n_m = (3^m n_0 + A) / 2^K.
    """
    if not word or any(int(k) < 1 for k in word):
        raise ValueError("word must be a nonempty sequence of positive valuations")
    m = len(word)
    Kprefix = 0
    A = 0
    for i, kval in enumerate(word):
        A += (3 ** (m - 1 - i)) * (1 << Kprefix)
        Kprefix += int(kval)
    return A, Kprefix


def valuation_word_residue(word: list[int] | tuple[int, ...]) -> tuple[int, int]:
    """Exact odd residue class realizing a finite accelerated valuation word.

    If A is the exact affine constant and K the total valuation, terminal
    oddness gives
        3^m n_0 + A == 2^K (mod 2^(K+1)).
    Since 3^m is invertible modulo 2^(K+1), this determines one odd
    residue class.  The classical parity-vector bijection guarantees that
    this congruence realizes every prescribed prefix valuation exactly.
    """
    A, K = valuation_word_affine_constant(word)
    modulus = 1 << (K + 1)
    inv = pow(pow(3, len(word), modulus), -1, modulus)
    residue = (((1 << K) - A) * inv) % modulus
    if residue == 0:
        residue = modulus
    return residue, modulus


def valuation_word_nondescend_threshold(word: list[int] | tuple[int, ...]) -> Fraction | None:
    """Exact threshold T=A/(2^K-3^m) for n_m >= n_0 when C>0.

    Returns None when 2^K <= 3^m, because terminal multiplicative
    compression is not positive.
    """
    A, K = valuation_word_affine_constant(word)
    den = (1 << K) - 3 ** len(word)
    if den <= 0:
        return None
    return Fraction(A, den)


def valuation_word_witness_count(word: list[int] | tuple[int, ...]) -> int | None:
    """Exact number of positive starts realizing word with C>0 but n_m >= n_0.

    The realizers are r+tM.  Intersecting this arithmetic progression with
    n <= A/(2^K-3^m) yields the closed-form finite count.  Returns None
    when C<=0 (the paradoxical/non-descent condition is not applicable).
    """
    threshold = valuation_word_nondescend_threshold(word)
    if threshold is None:
        return None
    residue, modulus = valuation_word_residue(word)
    if Fraction(residue, 1) > threshold:
        return 0
    return int((threshold - residue) // modulus) + 1


def valuation_word_witness_ratio(word: list[int] | tuple[int, ...]) -> float | None:
    """Return T/r; a positive-compression word has a witness iff T/r >= 1."""
    threshold = valuation_word_nondescend_threshold(word)
    if threshold is None:
        return None
    residue, _ = valuation_word_residue(word)
    return float(threshold / residue)

def realizes_valuation_word(n: int, word: list[int] | tuple[int, ...]) -> bool:
    """Return True iff odd n realizes word as its first accelerated valuations."""
    if n <= 0 or n % 2 == 0:
        return False
    x = int(n)
    for target in word:
        x, k = accelerated_step(x)
        if k != int(target):
            return False
    return True


def extremal_valuation_words(m: int, K: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return (back-loaded minimum-kernel word, front-loaded maximum-kernel word)."""
    if m < 1 or K < m:
        raise ValueError("require m >= 1 and K >= m")
    if m == 1:
        w = (K,)
        return w, w
    back = (1,) * (m - 1) + (K - m + 1,)
    front = (K - m + 1,) + (1,) * (m - 1)
    return back, front


def valuation_word_density_among_odds(word: list[int] | tuple[int, ...]) -> float:
    """Exact relative natural density within odd integers of a finite valuation word."""
    K = sum(int(k) for k in word)
    return 2.0 ** (-K)


def count_word_in_odd_progression(word, first_odd: int, nstarts: int) -> int:
    """Exact count of a word among first_odd, first_odd+2, ..., nstarts terms."""
    if nstarts < 0 or first_odd <= 0 or first_odd % 2 == 0:
        raise ValueError("first_odd must be positive odd and nstarts nonnegative")
    if nstarts == 0:
        return 0
    residue, modulus = valuation_word_residue(word)
    period = modulus >> 1  # index period after dividing the odd congruence by two
    target = ((residue - first_odd) // 2) % period
    if target >= nstarts:
        return 0
    return 1 + (nstarts - 1 - target) // period


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float); q = np.asarray(q, dtype=float)
    p = p / p.sum(); q = q / q.sum(); mid = 0.5 * (p + q)
    mp = p > 0; mq = q > 0
    return float(0.5 * np.sum(p[mp] * np.log(p[mp] / mid[mp])) +
                 0.5 * np.sum(q[mq] * np.log(q[mq] / mid[mq])))


def mutual_information_discrete(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x); y = np.asarray(y)
    ux, xi = np.unique(x, return_inverse=True)
    uy, yi = np.unique(y, return_inverse=True)
    table = np.zeros((len(ux), len(uy)), dtype=np.int64)
    np.add.at(table, (xi, yi), 1)
    pxy = table / table.sum(); px = pxy.sum(axis=1, keepdims=True); py = pxy.sum(axis=0, keepdims=True)
    den = px @ py; mask = pxy > 0
    return float(np.sum(pxy[mask] * np.log(pxy[mask] / den[mask])))


@dataclass
class BlockSummary:
    m: int
    mean_lambda: float
    sd_lambda: float
    anti_prob: float
    anti_prob_iid: float


def all_composition_count(m: int, K: int) -> int:
    """Number of positive valuation words of length m and total K."""
    if m < 1 or K < m:
        return 0
    return math.comb(K - 1, m - 1)


def witness_word_count_upper_bound(m: int, K: int) -> dict:
    """Rigorous residue-threshold upper bound for positive-compression witness words.

    A witness word is a positive valuation word of length m and total K for which
    at least one positive start realizes the word and nevertheless satisfies
    n_m >= n_0.  Distinct words have distinct odd representatives r modulo
    2^(K+1).  Since every witness obeys r <= T(word) <= T_max, the number of
    witness words is at most the number of positive odd integers <= T_max.

    This bound is exact as a counting argument and requires no stochastic
    independence assumption.
    """
    if m < 1 or K < m:
        raise ValueError("require m >= 1 and K >= m")
    den = (1 << K) - 3 ** m
    total = all_composition_count(m, K)
    if den <= 0:
        return {
            "m": m, "K": K, "compression_C": K * LN2 - m * LN3,
            "positive_compression": False, "T_max": math.inf,
            "witness_word_bound": total, "all_words": total,
            "fraction_bound": 1.0,
        }
    _, smax = residual_kernel_extrema(m, K)
    # A_max = 3^(m-1) S_max.  S_max is dyadic/rational but the float is
    # adequate for display; the integer-safe ceiling below is derived from
    # the exact front-loaded affine constant.
    front = (K - m + 1,) + (1,) * (m - 1)
    Amax, _ = valuation_word_affine_constant(front)
    # Number of positive odd integers r with r <= Amax/den is
    # floor((floor(Tmax)+1)/2).
    floor_t = Amax // den
    odd_count = (floor_t + 1) // 2
    bound = min(total, odd_count)
    return {
        "m": m, "K": K, "compression_C": K * LN2 - m * LN3,
        "positive_compression": True, "T_max": float(Fraction(Amax, den)),
        "witness_word_bound": int(bound), "all_words": int(total),
        "fraction_bound": float(bound / total),
    }


def strict_supercritical_witness_exponent(rho: float) -> dict:
    """Asymptotic exponents for K/m -> rho > log_2 3.

    The rigorous theorem in the manuscript gives
      log |W_{m,K}| <= m log(3/2) + o(m),
    while
      log C(K-1,m-1) = m*rho*H(1/rho)+o(m).
    Thus the witness-word fraction exponent is their difference.
    """
    critical = math.log(3, 2)
    if rho <= critical:
        raise ValueError("rho must be strictly greater than log_2 3")
    p = 1.0 / rho
    H = -p * math.log(p) - (1 - p) * math.log(1 - p)
    comp_exp = rho * H
    witness_exp = math.log(1.5)
    fraction_decay = comp_exp - witness_exp
    cylinder_decay = rho * LN2 - witness_exp
    return {
        "rho": rho, "composition_exponent": comp_exp,
        "witness_count_exponent_upper": witness_exp,
        "fraction_decay_exponent": fraction_decay,
        "odd_cylinder_decay_exponent": cylinder_decay,
    }


def compositions_positive(total: int, length: int):
    """Yield positive compositions; intended only for small exact checks."""
    if length == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - length + 2):
        for rest in compositions_positive(total - first, length - 1):
            yield (first,) + rest


def exact_witness_word_count_small(m: int, K: int, max_words: int = 2_000_000) -> dict:
    """Exhaustively count witness words for small (m,K), with a safety cap."""
    total = all_composition_count(m, K)
    if total > max_words:
        raise ValueError(f"composition family too large ({total:,} > {max_words:,})")
    if (1 << K) <= 3 ** m:
        raise ValueError("requires positive terminal compression 2^K > 3^m")
    witnesses = 0
    nearmax = 0
    max_ratio = 0.0
    for word in compositions_positive(K, m):
        ratio = valuation_word_witness_ratio(word)
        if ratio is not None:
            max_ratio = max(max_ratio, ratio)
            if ratio >= 1.0:
                witnesses += 1
    b = witness_word_count_upper_bound(m, K)
    return {
        "m": m, "K": K, "all_words": total,
        "exact_witness_words": witnesses,
        "exact_fraction": witnesses / total,
        "rigorous_bound": b["witness_word_bound"],
        "bound_fraction": b["fraction_bound"],
        "max_witness_ratio": max_ratio,
    }
