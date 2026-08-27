# Collatz Compression Spectra

## Reproducibility Software v2.0.0

**Copyright (C) 2026 Mohammad Amir Khusru Akhtar**  
Licensed under the [Apache License 2.0](LICENSE).

This repository provides the reproducibility software accompanying the manuscript:

**Compression Spectra of Collatz Orbits: Sharp Residual Ordering, Exact Valuation Cylinders, and Residue-Threshold Sparsity**

The software implements exact finite-time analysis of accelerated Collatz/Syracuse trajectories using Python arbitrary-precision integer arithmetic. It reproduces the computational tables and figures reported in the manuscript and provides an interactive application for exploring the theoretical framework under user-selected parameters.

## Main features

- exact accelerated Collatz/Syracuse iteration;
- arbitrary-precision integer arithmetic;
- compression and residual decomposition;
- fixed-`(m,K)` valuation-word analysis;
- sharp front-loaded and back-loaded residual extrema;
- universal non-descent, ordering-sensitive, and universal descent regimes;
- exact valuation-cylinder residue classes;
- residue-aware non-descent witness analysis;
- exact witness counts for tractable composition families;
- rigorous witness-word upper bounds;
- near-extremal valuation-word analysis;
- strict-supercritical asymptotic decay calculations;
- finite-block compression-spectrum experiments;
- excursion and stopping-time diagnostics;
- deterministic validation windows;
- manuscript-scale table and figure generation;
- automated theorem and consistency tests; and
- an interactive parameter-exploration application.

## One-click GitHub reproducibility

After uploading the contents of this folder to a GitHub repository:

1. Open **Actions**.
2. Select **One-click reproducibility**.
3. Click **Run workflow**.
4. Choose `paper` for the manuscript-scale analysis or `smoke` for a fast validation run.
5. Download the generated tables and figures from the workflow artifact.

The `paper_outputs/` directory contains the manuscript-scale tables and vector figures distributed with this release. These files are the same reproducibility artifact audited against the final manuscript.

## Interactive app

```bash
python -m pip install -r requirements.txt
python app.py
```

The app allows users to change the starting integer, accelerated block length `m`, total valuation `K`, explicit valuation word, near-extremal tolerance, and asymptotic ratio `rho`. It reports exact trajectory and compression/residual information, sharp extrema, start-size regimes, extremal residue classes, composition counts and densities, exact small-family witness counts where tractable, rigorous witness-word bounds, near-extremal bounds, and strict-supercritical decay exponents.

## Reproduce locally

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python run_analysis.py --out analysis_output --n-raw 200000 --n-spectrum 200000 --max-m 128 --n-validation 10000 --workers 2
```

## Scope

This software studies finite-time Collatz structure. Neither the software nor the accompanying manuscript claims a proof of the Collatz conjecture.

## Manuscript and citation

If you use this software or its accompanying methodology, please cite:

> Akhtar, M. A. K. (2026). *Compression Spectra of Collatz Orbits: Sharp Residual Ordering, Exact Valuation Cylinders, and Near-Extremal Sparsity* (Version V1). Zenodo. [https://doi.org/10.5281/zenodo.22015248](https://doi.org/10.5281/zenodo.22015248)

Machine-readable citation metadata: **[CITATION.cff](CITATION.cff)**

## License

This software is released under the **[Apache License 2.0](LICENSE)**.

Copyright (C) 2026 Mohammad Amir Khusru Akhtar.
