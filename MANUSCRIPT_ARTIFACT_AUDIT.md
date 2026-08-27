# Final Manuscript-Artifact Consistency Audit

This audit compares the final manuscript against the user-supplied `collatz-reproducibility-output.zip`. The manuscript tables are checked at their displayed precision; manuscript figures are copied byte-for-byte from the supplied artifact. Tiny last-bit floating-point differences seen in an earlier local software output do not affect any displayed manuscript value.

## Numerical tables and prose
- PASS - Table 1 m=8, K=13: artifact sharp_residual_extrema.csv
- PASS - Table 1 m=16, K=26: artifact sharp_residual_extrema.csv
- PASS - Table 1 m=32, K=51: artifact sharp_residual_extrema.csv
- PASS - Table 1 m=64, K=102: artifact sharp_residual_extrema.csv
- PASS - Table 2 m=8, K=13: artifact near_extremal_count_bounds.csv
- PASS - Table 2 m=16, K=26: artifact near_extremal_count_bounds.csv
- PASS - Table 2 m=32, K=51: artifact near_extremal_count_bounds.csv
- PASS - Table 2 m=64, K=102: artifact near_extremal_count_bounds.csv
- PASS - Table 3 m=16, K=28: artifact residue_threshold_word_bounds.csv
- PASS - Table 3 m=32, K=55: artifact residue_threshold_word_bounds.csv
- PASS - Table 3 m=64, K=109: artifact residue_threshold_word_bounds.csv
- PASS - Table 3 m=128, K=218: artifact residue_threshold_word_bounds.csv
- PASS - Table 4 (8, 13, 'front'): artifact count=24, normalized residue=0.66534423828125
- PASS - Table 4 (8, 13, 'back'): artifact count=25, normalized residue=0.01556396484375
- PASS - Table 4 (16, 26, 'front'): artifact count=0, normalized residue=0.9999949112534523
- PASS - Table 4 (16, 26, 'back'): artifact count=0, normalized residue=0.5942382737994194
- PASS - Table 4 (32, 51, 'front'): artifact count=0, normalized residue=0.6666666665890564
- PASS - Table 4 (32, 51, 'back'): artifact count=0, normalized residue=0.7503671646118162
- PASS - Table 4 (64, 102, 'front'): artifact count=0, normalized residue=1.0
- PASS - Table 4 (64, 102, 'back'): artifact count=0, normalized residue=0.24197510024714575
- PASS - Table 5 m=8: artifact compression_summary.csv
- PASS - Table 5 m=16: artifact compression_summary.csv
- PASS - Table 5 m=32: artifact compression_summary.csv
- PASS - Table 5 m=64: artifact compression_summary.csv
- PASS - Table 5 m=96: artifact compression_summary.csv
- PASS - Table 5 m=128: artifact compression_summary.csv
- PASS - Table 6 K=61: artifact exceptional_word_ordering.csv
- PASS - Table 6 K=65: artifact exceptional_word_ordering.csv
- PASS - Table 6 K=67: artifact exceptional_word_ordering.csv
- PASS - Table 6 K=70: artifact exceptional_word_ordering.csv
- PASS - Table 7 (17, 27): artifact accelerated_paradoxical_pairs.csv
- PASS - Table 7 (29, 46): artifact accelerated_paradoxical_pairs.csv
- PASS - Table 7 (41, 65): artifact accelerated_paradoxical_pairs.csv
- PASS - Table 7 (46, 73): artifact accelerated_paradoxical_pairs.csv
- PASS - Table 8 trajectory summary: artifact trajectory_summary.csv
- PASS - Table 9 holdout models: artifact multivariate_holdout_models.csv
- PASS - Table 10 cross-range: artifact out_of_sample_ranges.csv
- PASS - Table 11 generalized phase: artifact generalized_phase.csv
- PASS - Prose: m128 anti-count: matches artifact summary
- PASS - Prose: m128 exact density: matches artifact summary
- PASS - Prose: m128 tail rates: matches artifact summary
- PASS - Prose: deficit Spearman: matches artifact summary
- PASS - Prose: paradoxical census: matches artifact summary

## Figures
- PASS - `figures/residual_ordering_amplification.pdf` - SHA-256 `d6983b1fe4e5d46e380e74d17e393e0f03f7738b8e57a371c9cb7ef1d07c5ff0`
- PASS - `figures/residue_threshold_sparsity.pdf` - SHA-256 `5eb0a334b254ddaa492d6be038f275a55a24cebb9a664c8c967a1796c0534fb2`
- PASS - `figures/anti_compression_probability.pdf` - SHA-256 `1d03bbee086dd019422d7fbda37a5e3845203adaf7ed0c94de5f1dbbfa2acc73`
- PASS - `figures/tail_rate_zero.pdf` - SHA-256 `0c0d8ad6e457399a2f26a22b1d6f800060af07e76e5c14dc7b72ac33cb3d85b3`
- PASS - `figures/rate_function_collapse.pdf` - SHA-256 `c26ec99f20a90d1c882c320a961eb2bad609a4501f78e74ea30e229f3f96351d`
- PASS - `figures/compression_spectrum.pdf` - SHA-256 `ceb800faa8f933addc7ba0cbf78f3103abaf982afe4b7cc61844ea5dbf411712`
- PASS - `figures/valuation_distribution.pdf` - SHA-256 `0328103f028b619bed369a220b5bb721fbda2ff182791671ab5a56544ca70769`
- PASS - `figures/dependence_mi.pdf` - SHA-256 `292652d41ef3bd3bf3d216cac31cd9f41880066dbd27166afa04d5a44b60402a`
- PASS - `figures/deficit_vs_excursion.pdf` - SHA-256 `da5bf5b0e683a4fbfb86768c70705571838fc9031a539785ec233b8bd449cfcd`
- PASS - `figures/holdout_prediction.pdf` - SHA-256 `f8cc07ddda35d4153a98a7f79b94cffcb9e086e5eb0ca2702d440e58acb33451`
- PASS - `figures/cross_range_stability.pdf` - SHA-256 `e971cad99621121adcb6062476ca43ca017ab61d86d730a3922b787965fec215`
- PASS - `figures/stopping_survival.pdf` - SHA-256 `60ee62cd55c67545a60e89c5ce6251f9821a074c07109d56b3d9a7985a553b17`
- PASS - `figures/generalized_phase.pdf` - SHA-256 `66382eb0f8fd66e040ff42345f949125922a9a9ac21a68ceff2eee36991d9419`

## Citation coverage
- Bibliography entries: 70
- Unique cited entries: 70
- Uncited bibliography entries: none
- Missing bibliography keys: none

## Immediate float-reference audit
- PASS - table `tab:extrema` - immediately preceding nonblank line: `\Cref{tab:extrema} illustrates how rapidly the ordering-sensitive upper boundary can grow near the critical ratio $K/m\approx\log_2 3$. At $m=64$ and $K=102$, the same terminal compression permits a sharp residual-kernel ratio above $1.8\times10^{11}$, and the word-uniform descent threshold is approximately $3.84\times10^{11}$. Thus the terminal sum $K_m$ can be nearly uninformative about the worst possible residual until the starting value is sufficiently large.`
- PASS - figure `fig:orderingamp` - immediately preceding nonblank line: `\Cref{fig:orderingamp} visualizes the growth of the sharp ordering amplification for three near-critical $K$ choices at each scale.`
- PASS - table `tab:nearbound` - immediately preceding nonblank line: `\Cref{tab:nearbound} illustrates the quantitative bound at the same near-critical parameter choices used in \cref{tab:extrema}. Even at the relatively generous tolerance $\varepsilon=0.20$, the theorem allows at most one unit of excess valuation to be delayed for these cases, so the near-extremal family is at most linear in $m$ while the full composition class grows rapidly.`
- PASS - table `tab:witnessbound` - immediately preceding nonblank line: `\Cref{tab:witnessbound} reports the new arithmetic-filter bound at representative supercritical ratios. The bound is already tiny at moderate $m$ and decreases exponentially, as predicted by \cref{cor:witnessasymptotic}. \Cref{fig:witnesssparse} shows the same decay across several ratios.`
- PASS - figure `fig:witnesssparse` - immediately preceding nonblank line: `\Cref{fig:witnesssparse} visualizes the residue-threshold witness-word bound immediately following the numerical summary in \cref{tab:witnessbound}.`
- PASS - table `tab:residueclasses` - immediately preceding nonblank line: `\Cref{tab:residueclasses} makes the arithmetic sparsity concrete for the sharp extremizers used in \cref{tab:extrema}. At near-critical $K$, the density of either single extremal word decays as $2^{-K}$. Thus a finite interval may miss the front-loaded extremizer even though Corollary~\ref{cor:arithsharp} guarantees infinitely many global realizations.`
- PASS - algorithm `alg:certificate` - immediately preceding nonblank line: `Algorithms~\ref{alg:certificate} and \ref{alg:spectrum} summarize the two computational components of the framework before their detailed presentation: Algorithm~\ref{alg:certificate} provides the exact residue-aware certificate for an individual valuation word, whereas Algorithm~\ref{alg:spectrum} gives the reproducible multiscale compression-spectrum pipeline used to regenerate the empirical results.`
- PASS - algorithm `alg:spectrum` - immediately preceding nonblank line: `Algorithm~\ref{alg:spectrum} gives the manuscript-scale reproducibility pipeline used for all computational tables and figures.`
- PASS - table `tab:compression` - immediately preceding nonblank line: `\Cref{tab:compression} reports the multiscale results. Through $m=64$, empirical anti-compression frequencies lie close to the exact residue-density benchmark. At $m=128$, 16 events are observed among 200,000 starts, corresponding to $8.0\times10^{-5}$; the exact density is $8.862\times10^{-5}$, or 17.72 expected events. The Wilson 95\% interval for the empirical probability is $[4.92\times10^{-5},1.30\times10^{-4}]$, which contains the reference value.`
- PASS - figure `fig:anti` - immediately preceding nonblank line: `\Cref{fig:anti} compares the empirical and exact reference probabilities. The agreement remains close on a logarithmic scale, although the contiguous-window ensemble develops visible finite-window structure by $m=96$ and $128$.`
- PASS - figure `fig:tailrate` - immediately preceding nonblank line: `The tail-rate statistic $-m^{-1}\log P(\Lambda_m<0)$ is shown in \cref{fig:tailrate}. At $m=128$ the empirical value is $0.07370$ and the exact reference value is $0.07290$.`
- PASS - figure `fig:rate` - immediately preceding nonblank line: `for $m\in\{16,32,64,128\}$ where the empirical tail contains at least one observation. \Cref{fig:rate} compares these finite-scale curves with the Cram\'er reference \eqref{eq:cramer}. The curves retain lattice steps and window effects; therefore we do not claim an arithmetic large-deviation principle. The analysis does, however, turn ``compression spectrum'' into an explicit family of finite-scale rate curves with a falsifiable reference.`
- PASS - figure `fig:spectrum` - immediately preceding nonblank line: `The empirical distribution of $\Lambda_{128}$ is displayed in \cref{fig:spectrum}.`
- PASS - figure `fig:valuation` - immediately preceding nonblank line: `Pooling the complete 128-step valuation array gives Jensen--Shannon divergence $1.976\times10^{-3}$ nats from the clipped geometric reference. The larger divergence than in the earlier 32-step experiment is itself informative: long deterministic iteration of a narrow contiguous starting window need not behave like fresh independent residue sampling. \Cref{fig:valuation} shows the marginal comparison.`
- PASS - figure `fig:mi` - immediately preceding nonblank line: `Lagged mutual information is reported in \cref{fig:mi}. These diagnostics are deliberately descriptive: small pairwise information at short lags does not imply arithmetic independence, and departures at longer block scales can arise from residue geometry of the finite starting window.`
- PASS - table `tab:orderingreal` - immediately preceding nonblank line: `\Cref{tab:orderingreal} compares realized residual kernels at fixed $m=32$ and selected frequently occurring terminal valuations with the sharp theoretical all-word extrema from \cref{thm:ordering}. The realized words occupy only a tiny fraction of the worst-case kernel range. For example, among 10,296 realized words with $K=61$, the largest observed kernel is $6.28\times10^4$, whereas the all-word sharp maximum is $1.07\times10^9$.`
- PASS - table `tab:paradoxpairs` - immediately preceding nonblank line: `Proposition~\ref{prop:witnesscount} can be checked directly without probabilistic modeling. Scanning odd starts $3\le n_0\le200{,}000$ and accelerated prefixes through 100 odd events yields 148 positive-compression non-descent events from 144 distinct starts. Every event satisfies the exact residue congruence and threshold test in \eqref{eq:witnesscount}. \Cref{tab:paradoxpairs} reports the four $(m,K)$ pairs that occur. These compressed pairs correspond to odd-endpoint samples of the paradoxical finite behavior studied in parity-vector coordinates by \textcite{rozierterracol2026} and \textcite{niu2026}; the table is therefore a cross-formalism validation, not a new enumeration claim.`
- PASS - table `tab:trajectory` - immediately preceding nonblank line: `In the complete census $1\le n\le200{,}000$, every trajectory reached $1$ within the 10,000-step cap. The maximum total stopping time was 382 at $n=156159$, and the largest raw excursion was $17{,}202{,}377{,}752$ at $n=159487$ (\cref{tab:trajectory}).`
- PASS - figure `fig:deficit` - immediately preceding nonblank line: `The Spearman association between maximum compression deficit $D$ and $Y=\log(M(n)/n)$ is $\rho=0.669835$. Because the sample is very large, we emphasize effect size rather than a machine-underflowed $p$-value. \Cref{fig:deficit} displays the relationship.`
- PASS - table `tab:models` - immediately preceding nonblank line: `A deterministic random permutation with seed 20260819 assigns 70\% of the census to training and 30\% to holdout evaluation. Standardized OLS models are fitted only on the training partition. As \cref{tab:models} shows, $D$ alone achieves holdout $R^2=0.6360$. Adding the maximum residual changes this negligibly, while the full control model reaches $R^2=0.6465$. The standardized coefficient on $D$ remains $0.826$ in the full model, substantially larger than the coefficients on the other included variables.`
- PASS - figure `fig:holdout` - immediately preceding nonblank line: `\Cref{fig:holdout} compares observed and predicted holdout values for the full model.`
- PASS - table `tab:crossrange` - immediately preceding nonblank line: `The same analysis is repeated on three non-overlapping 10,000-start windows. \Cref{tab:crossrange} shows that the deficit--excursion rank association is stable between $0.677$ and $0.683$, while full-model holdout $R^2$ remains between $0.673$ and $0.707$. The standardized deficit coefficient remains between $0.864$ and $0.869$.`
- PASS - figure `fig:crossrange` - immediately preceding nonblank line: `\Cref{fig:crossrange} summarizes the stability visually.`
- PASS - figure `fig:survival` - immediately preceding nonblank line: `The stopping-time survival function in \cref{fig:survival} is shown descriptively; no heavy-tail claim is made without a separate asymptotic tail analysis \parencite{resnick2007}.`
- PASS - table `tab:phase` - immediately preceding nonblank line: `For odd $a$, multiplication by $a$ permutes odd residue classes modulo powers of two, so a fresh odd residue ensemble has $\E\nu_2(an+1)\approx2$. The multiplicative critical condition is therefore $\bar k>\log_2a$. \Cref{tab:phase,fig:phase} recover the familiar heuristic boundary between $a=3$ and $a=5$.`
- PASS - figure `fig:phase` - immediately preceding nonblank line: `\Cref{fig:phase} presents the same generalized-map calibration graphically.`

## Bibliographic corrections in this revision
- Added volume 99 to Alon, Behajaina, and Paran (2024), article 102473.
- Added volume 695 and pages 163-167 to Paparella (2024).
- Added and cited Paparella (2025) corrigendum, volume 708, pages 608-609, DOI 10.1016/j.laa.2024.12.019.

## Overall result
PASS if every item above is PASS. The supplied reproducibility artifact and manuscript agree at all reported/displayed precisions.
