import math
import itertools
import random
from collatz_compression import (
    accelerated_block, accelerated_step, iid_anticompression_probability,
    residual_kernel_from_word, residual_kernel_extrema, residual_bounds_fixed_total, near_extremal_count_bound,
    uniform_descent_certificate, valuation_word_residue, realizes_valuation_word, extremal_valuation_words, count_word_in_odd_progression, valuation_word_affine_constant, valuation_word_nondescend_threshold, valuation_word_witness_count, LN2, LN3
)


def compositions(total, parts):
    if parts == 1:
        yield (total,); return
    for first in range(1, total-parts+2):
        for rest in compositions(total-first, parts-1):
            yield (first,) + rest


def test_accelerated_step():
    assert accelerated_step(7) == (11, 1)
    assert accelerated_step(5) == (1, 4)


def test_exact_decomposition():
    n0 = 99991
    vals, ks, cs, rs = accelerated_block(n0, 20)
    for r in range(1, 21):
        lhs = math.log(vals[r] / vals[0])
        rhs = -cs[r-1] + rs[r-1]
        assert abs(lhs-rhs) < 2e-12


def test_iid_baseline():
    assert abs(iid_anticompression_probability(1)-0.5) < 1e-15
    assert abs(iid_anticompression_probability(4)-0.34375) < 1e-15


def test_residual_extrema_exhaustive_small():
    for m in range(2, 6):
        for K in range(m, m+6):
            vals = [residual_kernel_from_word(w) for w in compositions(K,m)]
            lo, hi = residual_kernel_extrema(m,K)
            assert abs(min(vals)-lo) < 1e-12
            assert abs(max(vals)-hi) < 1e-12



def test_near_extremal_count_bound_exhaustive_small():
    for m in range(2, 6):
        for K in range(m, m + 6):
            words=list(compositions(K,m))
            _,smax=residual_kernel_extrema(m,K)
            for eps in (0.0,0.05,0.20,0.30):
                info=near_extremal_count_bound(m,K,eps)
                actual=sum(residual_kernel_from_word(w) >= (1.0-eps)*smax-1e-12 for w in words)
                assert actual <= info['count_bound']
                assert info['count_bound'] <= len(words)
                if info['nontrivial']:
                    L=K-m
                    for w in words:
                        if residual_kernel_from_word(w) >= (1.0-eps)*smax-1e-12:
                            q=L-(w[0]-1)
                            assert q <= info['Q']

def test_uniform_certificate_consistency():
    n0=10**6+1; m=8; K=14
    lo,hi=residual_bounds_fixed_total(n0,m,K)
    C=K*LN2-m*LN3
    assert uniform_descent_certificate(n0,m,K) == (C>hi)


def test_exact_valuation_word_residues():
    for m in range(1, 5):
        for K in range(m, m + 5):
            seen=set()
            for w in compositions(K,m):
                residue,modulus=valuation_word_residue(w)
                assert residue % 2 == 1
                assert modulus == 2**(K+1)
                assert realizes_valuation_word(residue,w)
                assert residue not in seen
                seen.add(residue)
            assert len(seen) == math.comb(K-1,m-1)
            exact_density=len(seen) * 2.0**(-K)
            nb=math.comb(K-1,m-1)*2.0**(-K)
            assert abs(exact_density-nb) < 1e-15


def test_extremizers_are_arithmetically_realized():
    for m,K in [(2,4),(3,6),(4,7),(5,9)]:
        back,front=extremal_valuation_words(m,K)
        for w in {back,front}:
            residue,_=valuation_word_residue(w)
            assert realizes_valuation_word(residue,w)


def test_window_count_against_bruteforce():
    word=(1,1)
    first=10**4+1
    nstarts=100
    brute=sum(realizes_valuation_word(first+2*i,word) for i in range(nstarts))
    assert count_word_in_odd_progression(word,first,nstarts)==brute




def test_exact_witness_count():
    # The first nontrivial accelerated paradoxical prefix in the finite census.
    word=(4,1,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3)
    A,K=valuation_word_affine_constant(word)
    assert K==27
    residue,modulus=valuation_word_residue(word)
    assert residue==165
    assert realizes_valuation_word(165,word)
    T=valuation_word_nondescend_threshold(word)
    assert T is not None and float(T)>217 and float(T)<218
    assert valuation_word_witness_count(word)==1

if __name__ == '__main__':
    test_accelerated_step(); test_exact_decomposition(); test_iid_baseline()
    test_residual_extrema_exhaustive_small(); test_near_extremal_count_bound_exhaustive_small(); test_uniform_certificate_consistency(); test_exact_valuation_word_residues(); test_extremizers_are_arithmetically_realized(); test_window_count_against_bruteforce(); test_exact_witness_count()
    print('all core tests passed')


def test_witness_word_count_bound_small():
    from collatz_compression import exact_witness_word_count_small, witness_word_count_upper_bound
    for m, K in [(3,5),(4,7),(5,8),(6,10)]:
        if (1 << K) <= 3**m:
            continue
        ex = exact_witness_word_count_small(m,K,max_words=500000)
        bd = witness_word_count_upper_bound(m,K)
        assert ex['exact_witness_words'] <= bd['witness_word_bound']


def test_strict_supercritical_exponent_positive():
    from collatz_compression import strict_supercritical_witness_exponent
    out = strict_supercritical_witness_exponent(1.7)
    assert out['fraction_decay_exponent'] > 0
    assert out['odd_cylinder_decay_exponent'] > 0
