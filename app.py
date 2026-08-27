"""Interactive Gradio explorer for the Collatz compression framework."""
from __future__ import annotations
import math
import pandas as pd
import gradio as gr
from collatz_compression import (
    accelerated_block, residual_kernel_extrema, residual_bounds_fixed_total,
    extremal_valuation_words, valuation_word_residue,
    valuation_word_nondescend_threshold, valuation_word_witness_count,
    valuation_word_witness_ratio, near_extremal_count_bound,
    witness_word_count_upper_bound, strict_supercritical_witness_exponent,
    iid_anticompression_probability, all_composition_count,
    exact_witness_word_count_small, LN2, LN3,
)


def trajectory_view(n, m):
    n, m = int(n), int(m)
    vals, ks, cs, rs = accelerated_block(n, m)
    df = pd.DataFrame({
        "step": range(1, m + 1), "odd_state": vals[1:], "valuation_k": ks,
        "compression_C": cs, "residual_R": rs, "log_displacement": -cs + rs,
    })
    summary = {
        "terminal_compression": float(cs[-1]), "terminal_residual": float(rs[-1]),
        "terminal_log_displacement": float(-cs[-1] + rs[-1]),
        "anti_compression": bool(cs[-1] < 0),
    }
    return df, summary


def fixed_mk_view(m, K, n0, epsilon):
    m, K, n0 = int(m), int(K), int(n0)
    if K < m:
        return {"error": "K must be at least m"}, pd.DataFrame()
    smin, smax = residual_kernel_extrema(m, K)
    rmin, rmax = residual_bounds_fixed_total(n0, m, K)
    C = K * LN2 - m * LN3
    back, front = extremal_valuation_words(m, K)
    total_words = all_composition_count(m, K)
    out = {
        "m": m, "K": K, "K_over_m": K / m,
        "C": C, "terminal_multiplicative_factor_3m_over_2K": math.exp(-C),
        "all_positive_valuation_words": total_words,
        "exact_relative_density_K_among_odd_starts": total_words * (2.0 ** (-K)),
        "S_min": smin, "S_max": smax, "S_max/S_min": smax/smin,
        "R_min_at_n0": rmin, "R_max_at_n0": rmax,
        "back_loaded_word": back, "front_loaded_word": front,
        "back_loaded_residue": valuation_word_residue(back),
        "front_loaded_residue": valuation_word_residue(front),
        "anti_compression_density_at_m": iid_anticompression_probability(m),
    }
    if C > 0:
        den = 3.0 * math.expm1(C) if C < 700 else math.inf
        n_minus = smin / den if math.isfinite(den) else 0.0
        n_plus = smax / den if math.isfinite(den) else 0.0
        out.update({
            "n_minus_universal_nondescend_boundary": n_minus,
            "n_plus_universal_descent_boundary": n_plus,
            "n0_regime": ("universal non-descent" if n0 < n_minus else
                          "universal descent" if n0 > n_plus else
                          "ordering-sensitive strip"),
        })
    else:
        out["n0_regime"] = "no positive terminal compression"

    witness = witness_word_count_upper_bound(m, K)
    rows = [{"result": "witness-word rigorous upper bound", **witness}]
    if total_words <= 200_000 and C > 0:
        exact = exact_witness_word_count_small(m, K, max_words=200_000)
        rows.append({"result": "exact witness-word census (tractable case)", **exact})
    if m >= 2:
        near = near_extremal_count_bound(m, K, float(epsilon))
        rows.append({"result": "near-extremal count bound", **near})
    return out, pd.DataFrame(rows)


def word_view(raw):
    try:
        word = tuple(int(x.strip()) for x in raw.split(",") if x.strip())
        if not word or any(k < 1 for k in word):
            raise ValueError
        residue, modulus = valuation_word_residue(word)
        T = valuation_word_nondescend_threshold(word)
        out = {
            "word": word, "m": len(word), "K": sum(word), "residue": residue,
            "modulus": modulus, "conditional_density_among_odds": 2.0 ** (-sum(word)),
        }
        if T is None:
            out["status"] = "No positive terminal multiplicative compression (2^K <= 3^m)."
        else:
            out.update({"threshold_T": str(T), "witness_ratio_T_over_r": valuation_word_witness_ratio(word),
                        "exact_nondescend_start_count": valuation_word_witness_count(word)})
        return out
    except Exception as exc:
        return {"error": f"Invalid valuation word: {exc}"}


def asymptotic_view(rho):
    try:
        return strict_supercritical_witness_exponent(float(rho))
    except Exception as exc:
        return {"error": str(exc)}


with gr.Blocks(title="Collatz Compression Spectra Explorer") as demo:
    gr.Markdown("# Collatz Compression Spectra Explorer\nExact finite-time theorem and diagnostic explorer. **No output is a proof of the Collatz conjecture.**")
    with gr.Tab("Start trajectory"):
        with gr.Row():
            n = gr.Number(value=27, label="Positive starting integer", precision=0)
            m = gr.Slider(1, 256, value=32, step=1, label="Accelerated odd steps")
        btn = gr.Button("Compute exact trajectory")
        table = gr.Dataframe(label="Exact path data")
        summary = gr.JSON(label="Summary")
        btn.click(trajectory_view, [n, m], [table, summary])
    with gr.Tab("Fixed (m,K)"):
        with gr.Row():
            m2 = gr.Number(value=16, label="m", precision=0)
            K2 = gr.Number(value=28, label="K", precision=0)
            n0 = gr.Number(value=1000001, label="Reference n0", precision=0)
            eps = gr.Slider(0, .95, value=.20, step=.01, label="Near-extremal epsilon")
        btn2 = gr.Button("Compute fixed-(m,K) bounds")
        summary2 = gr.JSON(label="Sharp bounds")
        table2 = gr.Dataframe(label="Residue-threshold / near-extremal bounds")
        btn2.click(fixed_mk_view, [m2,K2,n0,eps], [summary2,table2])
    with gr.Tab("Valuation word"):
        raw = gr.Textbox(value="4,1,1,1,1", label="Comma-separated positive valuations")
        btn3 = gr.Button("Certify word")
        out3 = gr.JSON(label="Exact residue-aware certificate")
        btn3.click(word_view, raw, out3)
    with gr.Tab("Asymptotic theorem"):
        rho = gr.Number(value=1.70, label="rho (> log2 3)")
        btn4 = gr.Button("Compute decay exponents")
        out4 = gr.JSON(label="Strict-supercritical exponents")
        btn4.click(asymptotic_view, rho, out4)

if __name__ == "__main__":
    demo.launch()
