# WP2(a) — alpha-claim audit: which manuscript statements about the MC-SRT significance level survive the NN-table duplication

Driver: `revision_experiments/66_alpha_claim_audit.R`
Per-claim CSV: `revision_experiments/results/tr1/wp2a_alpha_claim_audit.csv` (24 rows: 16 real data sets + 7 simulation dimensions + header)

This is a purely documentary audit — no Monte Carlo was run. It takes as given
`62_audit_nn_duplicates.R`'s finding that 21 (d, token-pair) combinations in
`R/NN-test_quantile/` are byte-identical files saved under two alpha labels
(d = 6,7,8,9 for the 95%/99% pair; d = 11–19, 21–28 for the 99%/99.9% pair),
and asks a narrower, mechanical question: given the code's *actual* token
resolvers, which of the manuscript's specific alpha-level sentences rest on a
distinct file, and which rest on a duplicate?

## Headline

- **No claim in the manuscript is FALSE in the "code loads a different label
  than stated" sense (mismatch type (c)).** Across all 48 real-data claims
  (16 data sets × RK/UN/SUN) and all 21 simulation-grid claims (7 dimensions
  × RK/UN/SUN), the token the code resolves always matches the percentage the
  manuscript states. Row-range placement (type (d)) is also clean: all 16 data
  sets sit in the table row their actual `d` puts them in.
- **The RK column of Table `tab:alpha_real` is fully sound.** Every RK file
  used by the 16 real data sets is either the only token that exists at that
  `d` (12, 16, 18, 19, 21, 30) or is byte-provably distinct from its sibling
  at the same `d` (5, 6, 7, 8, 9, 10, 20 — file sizes differ by 66 KB to 326
  MB, which alone rules out byte-identity; no RK pair anywhere in the
  dimensions used by this paper needed a content load to settle).
- **The simulation schedule (SM L718–740, main text L907–910) is fully
  sound**, including the verification comment at SM L740. All 21 claims at
  `d ∈ {2,3,5,10,20,50,100}` are VERIFIED; the d = 10 UN-MCCD/SUN-MCCD split
  the comment asserts is real (`NN-test-simul_10d_99%.RData` is 70717 bytes,
  `NN-test-simul_10d_999%.RData` is 70940 bytes — different sizes, and a full
  `$average`/`$median` load confirms different content).
- **The real-data NN column (UN-MCCD, SUN-MCCD) is where the damage is.** Of
  the 32 real-data NN claims (16 data sets × 2 methods), **26 are
  UNVERIFIABLE** — the loaded file is a content-duplicate of the sibling
  token at that `d`, so its true generating alpha cannot be distinguished
  from the sibling's. Only wilt (d=5), pageblocks (d=10), and WDBC (d=30,
  vacuously — no sibling file exists to duplicate) are clean.
- **Four data sets carry a claim that is affirmatively FALSE, not merely
  unverifiable**: hepatitis (d=19), vowels (d=12), lymphography (d=18), and
  PenDigits (d=16). Table `tab:alpha_real`'s row 2 states UN-MCCD uses 1% and
  SUN-MCCD uses 0.1% for these — different numbers — but the code loads
  **byte-identical files** for both methods at each of these four dimensions.
  The claimed distinction did not exist in the computation that produced the
  paper's numbers.

## Per-data-set verdict table

RK is sound throughout, so only the NN (UN-MCCD/SUN-MCCD) verdicts vary.
"file(s)" gives the loaded NN filename(s); "α claimed" is Table
`tab:alpha_real`'s value; "duplicate of" names the sibling token the file is
byte-identical to, where applicable.

| Data set | n | d | row claimed | UN α claimed | UN verdict | SUN α claimed | SUN verdict | UN vs SUN distinction |
|---|---|---|---|---|---|---|---|---|
| wilt | 4819 | 5 | 5–9 | 5% | **VERIFIED** — distinct from the 90% file | 5% | **VERIFIED** | same token, none claimed |
| vertebral | 240 | 6 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| thyroid | 3656 | 6 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| ecoli | 336 | 7 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| pima | 555 | 8 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| glass | 213 | 9 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| WBC | 223 | 9 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| Shuttle | 1013 | 9 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| stamps | 340 | 9 | 5–9 | 5% | UNVERIFIABLE — dup. of 99%(1%) | 5% | UNVERIFIABLE | same token, none claimed |
| pageblocks | 4795 | 10 | 10–19 | 1% | **VERIFIED** — distinct from 99.9% file | 0.1% | **VERIFIED** — distinct | **real: files genuinely differ (70717 vs 70940 bytes)** |
| vowels | 1452 | 12 | 10–19 | 1% | UNVERIFIABLE — dup. of 999%(0.1%) | 0.1% | UNVERIFIABLE — dup. of 99%(1%) | **FALSE — same file loaded for both, claimed levels differ** |
| PenDigits | 3200 | 16 | 10–19 | 1% | UNVERIFIABLE — dup. of 999%(0.1%) | 0.1% | UNVERIFIABLE — dup. of 99%(1%) | **FALSE — same file loaded for both** |
| lymphography | 148 | 18 | 10–19 | 1% | UNVERIFIABLE — dup. of 999%(0.1%) | 0.1% | UNVERIFIABLE — dup. of 99%(1%) | **FALSE — same file loaded for both** |
| hepatitis | 74 | 19 | 10–19 | 1% | UNVERIFIABLE — dup. of 999%(0.1%) | 0.1% | UNVERIFIABLE — dup. of 99%(1%) | **FALSE — same file loaded for both** |
| waveform | 3443 | 21 | ≥20 | 0.1% | UNVERIFIABLE — dup. of 99%(1%) | 0.1% | UNVERIFIABLE — dup. of 99%(1%) | same token, none claimed |
| WDBC | 367 | 30 | ≥20 | 0.1% | VERIFIED* — single file, no sibling to cross-check | 0.1% | VERIFIED* | same token, none claimed |

\* WDBC's verdict is the weakest form of "VERIFIED": only one NN token file
exists at d=30, so there is nothing to compare it against and no way to
detect a mislabeling the way the other 21 duplicate pairs were caught. Treat
it as "no contrary evidence," not as independently confirmed.

## Three-way breakdown (48 real-data claims: 16 data sets × {RK, UN, SUN})

| | VERIFIED | UNVERIFIABLE | FALSE (type c: wrong token) |
|---|---|---|---|
| RK (16 claims) | 16 | 0 | 0 |
| UN (16 claims) | 3 (wilt, pageblocks, WDBC) | 13 | 0 |
| SUN (16 claims) | 3 (wilt, pageblocks, WDBC) | 13 | 0 |
| **Total** | **22** | **26** | **0** |

Plus, orthogonal to the per-cell table above, **4 of the 5 data sets in the
"10–19" row (all but pageblocks) carry a FALSE UN-vs-SUN *distinction*
claim** (case (b) in the task brief): the table states different numbers for
UN-MCCD and SUN-MCCD, but the code loads one file for both.

Simulation grid (21 claims: 7 dimensions × {RK, UN, SUN}): **21 VERIFIED, 0
UNVERIFIABLE, 0 FALSE.** The duplication problem does not touch the
simulation schedule.

## What "VERIFIED" and "UNVERIFIABLE" mean here, precisely

- **VERIFIED**: the token the code resolves matches the manuscript's stated
  percentage, *and* the resulting file's content (`$average`/`$median` for
  NN, `$quan` for RK) is provably different from every other token's file at
  the same `d` — so the file is not a mislabeled copy of a table meant for a
  different alpha.
- **UNVERIFIABLE**: token matches, file exists, but it is byte-identical to
  a file saved under a *different* alpha label at the same `d`. One
  Monte-Carlo run produced both files; only one of the two labels can be
  correct, and there is no remaining evidence in the repository that says
  which. The manuscript's stated percentage for these cells is therefore
  neither confirmed nor contradicted by the file itself — it could be right,
  but nothing on disk proves it.
- **FALSE** (reserved for case (c): code resolves a token whose percentage
  disagrees with what the manuscript states) — **did not occur anywhere in
  this audit.** The separate FALSE finding that did occur is the
  UN-vs-SUN-distinction claim (case (b)), which is a claim about a
  *relationship between two cells*, not about a single cell's label.

## Manuscript locations requiring correction

### 1. `CCD_OutlierDetection_Neurocomputing.tex`, Table `tab:alpha_real`, lines 1104–1110

Current text:

```
1104:  \textbf{Data sets} & $d$ & \textbf{U-MCCD, SU-MCCD} & \textbf{UN-MCCD} & \textbf{SUN-MCCD} \\ \hline
1105:  wilt, vertebral, thyroid, ecoli, pima, glass, WBC, Shuttle, stamps & $5$--$9$   & $1\%$   & $5\%$   & $5\%$   \\ \hline
1106:  pageblocks, vowels, PenDigits, lymphography, hepatitis             & $10$--$19$ & $0.1\%$ & $1\%$   & $0.1\%$ \\ \hline
1107:  waveform, WDBC                                                     & $\geq 20$  & $0.1\%$ & $0.1\%$ & $0.1\%$ \\ \hline
1108:\end{tabular}}
1109: \caption{MC-SRT significance level $\alpha$ used for each real data set, obtained by applying the dimension-specific schedules of Section~\ref{sec:Simul_CCDs} as a step function. The RK-based methods switch at $d=10$; UN-MCCD switches at $d=10$ and again at $d=20$; SUN-MCCD switches once, at $d=10$.}
```

Problem: the UN-MCCD and SUN-MCCD columns (not U-MCCD/SU-MCCD, which are
sound) present single numbers as established facts for every listed data
set. For 13 of the 16 data sets that number is unverifiable (row-1 data sets
except wilt, plus waveform), and for 4 of the 16 (vowels, PenDigits,
lymphography, hepatitis) the table additionally implies UN-MCCD and SUN-MCCD
ran at different significance levels when, in the actual computation, they
consulted the same file. Only wilt, pageblocks, and WDBC in the UN/SUN
columns are on solid footing (and WDBC only because no rival file exists to
contradict it).

Needed correction (a decision for the user/paper, not made here): either (a)
add a footnote to the table disclosing that the NN quantile tables for the
listed dimensions carry duplicate-content pairs, so the stated NN α values at
those dimensions could not be independently distinguished from their
step-function neighbor and the vowels/PenDigits/lymphography/hepatitis
UN-vs-SUN split did not exist computationally, similar in spirit to the
`\revised{}` disclosure already given to the withdrawn saturation claim; or
(b) regenerate the affected NN tables with a real Monte Carlo run so the
claimed values become verifiable, then update the table. Given the CLAUDE.md
protection rule (numerical results are read-only without explicit sign-off),
this audit stops at identifying the problem.

### 2. `CCD_OutlierDetection_Neurocomputing.tex`, lines 907–912 (prose)

Current text:

```
907: For the RK-based methods, U-MCCD and SU-MCCD, we use $\alpha=1\%$ for $d<10$ and $\alpha=0.1\%$ for the MC-SRT when $d\geq 10$. 
908: For the NND-based methods, UN-MCCD and SUN-MCCD, the significance level $\alpha$ is dimension-specific.
909: UN-MCCD uses $\alpha=15\%,10\%,5\%,1\%,0.1\%$ for $d=2,3,5,10,\{20,50,100\}$, respectively,
910: while SUN-MCCD uses $\alpha=0.1\%$ already at $d=10$; see the Supplementary Material for details.
911: These levels are specified at the tabulated dimensions and are applied as a step function elsewhere: each level holds for every dimension from its own anchor up to, but not including, the next one.
912: The rule matters for the real data of Section~\ref{sec:Real-Data-Examples}, whose dimensions do not coincide with the simulation grid; Table~\ref{tab:alpha_real} lists the resulting levels.
```

This prose describes the *simulation-grid* schedule, which is fully verified
(see below) — no correction needed to L907–910 themselves. L912 is the
sentence that hands the reader off to Table `tab:alpha_real` and calls its
contents "the resulting levels"; once the table above carries a duplicate-data
caveat, this sentence should point to it (or the caveat should be attached
here instead of at the table — a phrasing choice, not something this audit
should decide).

### 3. `SupplementaryMaterial.tex`, lines 717–740 — no correction needed

Quoted for completeness, since it is the schedule the real-data table claims
to derive from:

```
718: The MC-SRT significance level is set to $\alpha=1\%$ for $d<10$ and $\alpha=0.1\%$ for $d\geq 10$. 
719: For the NND-based methods, the dimension-specific levels coincide for $d\in\{2,3,5\}$ and for $d\in\{20,50,100\}$, and differ only at $d=10$. For UN-MCCD,
...
740: % Verified against the released simulation drivers: each per-dimension script loads a precomputed NND quantile file (R/NN-test_quantile/) that fixes the test level. At d=10, UN-MCCD loads the 99% file (alpha=1%) and SUN-MCCD loads the 99.9% file (alpha=0.1%); all other dimensions coincide.
```

**This assertion still holds.** All 7 simulation dimensions × 3 methods (21
claims) are VERIFIED by this audit, including the specific d=10 UN/SUN split
the comment calls out: `NN-test-simul_10d_99%.RData` (70,717 bytes) and
`NN-test-simul_10d_999%.RData` (70,940 bytes) have different sizes and,
confirmed by a full load, different `$average`/`$median` vectors. Do not
touch this file on the strength of this audit — the problem is confined to
how the *real-data* dimensions (which don't coincide with the simulation
anchors) land on the NN table directory, not to the simulation schedule
itself.

## What this audit could not determine

- **Which of the two labels is correct for each duplicated file.** E.g. for
  `NN-test-simul_9d_95%.RData` / `NN-test-simul_9d_99%.RData` (byte-identical),
  nothing in the repository says whether the underlying Monte Carlo run was
  actually configured for α=5% or α=1% — both filenames were applied to one
  run. `63_nn_provenance_recon.R` may have relevant provenance evidence but
  reading it was out of scope for this mechanical audit.
- **Whether the outlier flags themselves change** if the true alpha differs
  from the claimed one. This audit only established which *files* are
  distinguishable, not how sensitive each data set's BA/F2 would be to
  loading the "other" file — that is the subject of the pre-existing
  `wp2a_alpha_sweep.R` / `wp2a_alpha_findings.md` work, which found (for the
  RK side, where regeneration is free) that BA moves by a median of 0.16
  across the tested alpha ladder and is non-monotone in 19 of 22 groups.
  Because the NN side cannot be regenerated without a new Monte Carlo run,
  no comparable sensitivity number exists for UN-MCCD/SUN-MCCD.
- **WDBC's single-file dimension (d=30)** cannot be cross-checked at all —
  there is no sibling token file to compare against, so a mislabeling there
  (if one exists) is undetectable by this method.
