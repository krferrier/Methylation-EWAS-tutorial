# Resolved decisions for the v8.1 rebuild

Settled before any recompute. Each entry records what was decided, by whom, and the
evidence behind it.

## 1. Mask source: Zhou v8.1 `EPIC.hg38.mask.cm`, read with YAME

**Decided by:** user.
**Instruction:** use the most recent Zhou EPIC `mask.cm`, not the legacy
`EPIC.hg38.manifest.pop.tsv.gz`; note that the legacy file with sub-population
masking remains available where it is more appropriate.

The `.cm` file is not in `zhou-lab/InfiniumAnnotationData` (which I searched first,
and which carries only `EPIC.hg38.mask.tsv.gz`). It lives in
`zhou-lab/InfiniumAnnotation` at `EPIC/EPIC.hg38.mask.cm`, and is fetched by YAME
itself. The YAME catalog labels that directory **v8.1**, confirming the version
named in the revised prose.

The `.cm` is a bit-packed CX-format file: 24 mask terms over the full
**866,553-probe** EPIC ordering, one bit per probe, 357 KB for mask + index. This is
a materially better teaching object than the TSV, which lists only already-masked
probes and requires parsing a comma-delimited `maskUniq` column.

**YAME provenance:** built from source at **v1.40**. Bioconda's `linux-64` build is
pinned at 1.8, which predates `yame fetch` and so cannot download the reference data
that `-R`/`-m` resolve against. The `zhou-lab` conda channel (v1.40) is not on this
environment's channel allowlist, so source build was the route. Archived as
`yame_v140_build_and_store.tar.gz`.

## 2. Mask columns: the eight named codes

**Decided by:** user.

`M_mapping`, `M_nonuniq`, `M_SNP_AFR_1pt`, `M_1baseSwitchSNP_AFR_1pt`,
`M_2extBase_SNP_AFR_1pt`, `M_SNP_EUR_1pt`, `M_1baseSwitchSNP_EUR_1pt`,
`M_2extBase_SNP_EUR_1pt` — AFR and EUR both included because the Grady cohort is
~95% African American and admixed.

Per-term counts read from the `.cm` (`yame unpack -a`, cross-checked against
`yame summary`):

| term | probes |
|---|---|
| `M_mapping` | 428 |
| `M_nonuniq` | 28,833 |
| `M_SNP_AFR_1pt` | 35,420 |
| `M_1baseSwitchSNP_AFR_1pt` | 445 |
| `M_2extBase_SNP_AFR_1pt` | 20,640 |
| `M_SNP_EUR_1pt` | 17,076 |
| `M_1baseSwitchSNP_EUR_1pt` | 259 |
| `M_2extBase_SNP_EUR_1pt` | 8,887 |
| **union** | **90,832** |

**Validation:** the union was computed independently twice — from the bit-packed
`.cm` via YAME, and from `EPIC.hg38.mask.tsv.gz` via `data.table` — and the two probe
sets are `identical()`. All 90,832 fall inside the tutorial's EPIC v1 probe universe
(0 outside).

For reference, the as-run (superseded) filtering used `MASK_general_AFR` from the
legacy pop manifest: 124,132 probes, retaining 724,282.

## 3. Recompute scope: chapters 03 through 08

**Decided by:** user — "a full recompute of 03 through 08 should be done for complete
internal consistency."

**Chapter 04 is excluded, deliberately.** `estimateCellCounts2` runs on the raw
`RGChannelSet` with `processMethod="preprocessNoob"` and IDOL-optimized probes, so
cell composition is computed upstream of probe filtering and is unaffected by it.
Its checkpoints stay as-is.

## 4. Probe scope: all probes passing QC and filtering

**Decided by:** user — "I would like the tutorial to test all probes that pass QC and
filtering, not just the 20,000."

Two distinct 20,000-probe caps existed, and neither was the limma EWAS (which was
already full-scale at 724,282 probes):

- **ch05** — PCA computed on the 20,000 most variable probes (`fig-cap`, line 58).
- **ch07** — a genuine `set.seed(42); sample(rownames(betas), 20000)` at line 109,
  with prose and two runtime claims ("~11 s", "~15 seconds") that will not survive
  the change of scale.

Both are removed. Runtimes at full scale will be measured, not predicted.

## 5. Cell composition: keep all six proportions

**Decided by:** agent, by direct test, rather than referred to the user.

The revised 06 prose uses all 6 proportions where the previous version dropped one to
avoid compositional collinearity, and deleted the paragraph explaining why. I tested
whether that matters:

| | 6 proportions | 5 proportions |
|---|---|---|
| design columns | 25 | 24 |
| rank | 25 | 24 |
| full rank | TRUE | TRUE |
| residual df | 62 | 63 |
| VIF on PTSD coefficient | 1.224 | 1.224 |
| condition number | 25,400 | 7,486 |

The proportions sum to 0.982–1.048 in this data (only 1 of 96 samples sums to exactly
1), so they are not exactly compositional and the design does not become singular.

**Conclusion:** the 6-proportion design is estimable. It costs one degree of freedom
and raises the condition number, but leaves the precision of the PTSD coefficient
unchanged — not the harm the deleted paragraph implied. Keeping 6 as the prose states,
and describing the trade-off accurately rather than restoring the old claim.

## 6. Smoking

**Decided by:** user. Smoking is a major methylation variance component; the Grady
deposit provides no smoking variable, and many public datasets will not either. Since
SVA is already used and no smoking term is in the null model, smoking-driven variation
is expected to be captured by the surrogate variables, and that is how this tutorial
addresses it. The fuller option list belongs in 06; preprocessing chapters get a
cross-referencing mention.

## 7. CGI figure

**Decided by:** user — supplying `data/CGI.png` directly; no schematic to be drawn.
The surrounding prose needs its region distances reconciled against the supplied
image (shore 0–2 kb, shelf 2–4 kb, but open sea currently stated as "greater than
2kb", which overlaps the shelf).

## Smoking proxy, SVA null model, and ComBat (user rulings, this session)

### The proxy goes in the SVA null model, not merely alongside the SVs
User: "Couldn't the smoking proxy be included in the null model for SVA or used
to protect against with ComBat? Also, any CpGs used for making the smoking proxy
should be excluded from the EWAS."

Rationale: SVA estimates surrogate variables from variation *not* explained by
the supplied model. Placing the proxy in both `mod` and `mod0` protects it, so
SVA no longer reconstructs smoking as a surrogate variable and the SVs represent
only residual unmeasured variation. Smoking becomes an explicitly modeled
covariate rather than something hoped to be captured.

### Panel CpGs are excluded from the tested set
The 20 panel probes present after v8.1 filtering are removed from the M-value
matrix before SVA and before the EWAS. A probe cannot both define a covariate
and be tested against it. This matches the option already recommended in the
chapter 01 prose ("prior is recommended if using this option"). Cost: 20 probes.

### Smoking proxy construction (for the tutorial's methods text)
1. Panel of 24 replicated smoking CpGs across AHRR, F2RL3, GPR15, 2q37.1,
   GFI1, MYO1G, PRSS23, RARA.
2. 20 of 24 survive v8.1 masking; absent: cg23576855, cg05951221, cg06126421,
   cg22132788. The anchor cg05575921 survives.
3. Score = PC1 of the row-z-scored panel submatrix (52.9% of panel variance).
4. Sign oriented on AHRR cg05575921 (hypomethylated in smokers); anchor
   correlation -0.908. All loadings negative except the three MYO1G probes.
Refs: joehanes2016smoking, elliott2014smoking, breitling2011f2rl3,
bauer2015gpr15, shenker2013epic, zeilinger2013smoking. Alternatives to name:
mccartney2018scores, zhang2016mortality (weighted scores),
bollepalli2019epismoker (trained classifier).

### Why the proxy does NOT replace the SVs
R2(proxy ~ 15 SVs) = 0.625 under the pre-proxy model: the SVs captured ~63% of
the proxy, leaving 37% unexplained -- neither substitutes for the other. The SVs
also carry technical structure the proxy is blind to (SV4 chip R2 = 0.88).

### ComBat: batch structure of GSE132203 subset
R2(PTSD ~ slide)    = 0.058, chisq p = 0.95   -> PTSD is balanced across slides
R2(sex  ~ slide)    = 0.874, chisq p = 4.9e-5 -> slide nearly nested within sex
R2(PTSD ~ position) = 0.075
5 slides are all-male, 5 all-female. ComBat on slide is therefore admissible
for PTSD but MUST protect sex in its covariate model, or it removes real sex
signal. See nygaard2016batch for this failure mode and for the anticonservative
p-values that follow when ComBat-corrected data is passed downstream as if
untouched -- relevant at n=87 with a 25-parameter design.

### comb-p DMR results now exist (was blocked)
Fixed by shadowing `toolshed` with a patched copy: `threading.TIMEOUT_MAX` does
not exist in Python 2.7 (replaced with 2147483647) -- see `pyshim/`.
Strict (--seed 1e-15): 1 region, chr6:28,633,587-28,633,601, 3 probes,
  z_sidak_p = 9.774e-06. ACF corr 0.0527 out to 71 bp over 202,084 pairs.
Demo (--seed 1e-8): 1 region, chr6:28,633,493-28,633,667, 10 probes,
  z_sidak_p = 1.741e-19.
NOTE: `min_p` in the regions output is the SLK-smoothed p, NOT a raw probe p --
the strict region reports 1.516e-17 while the smallest raw BACON p in the BED is
1.63e-08. The prose must not present min_p as a probe-level p-value.
Final annotation step needs `bedtools` (was absent; installed into env combp).

### Upstream bug found: missMethyl 1.32.0 KEGG species code
`missMethyl:::.getKEGG()` calls `getGeneKEGGLinks(species.KEGG = "hsa")` but
`getKEGGPathwayNames(species.KEGG = "Hsa")`. KEGG's REST API is case-sensitive
and returns HTTP 400 for "Hsa", so `gometh(collection = "KEGG")` fails for every
user of this version. Patched in-session in recompute_08_enrich.R.

## Stratified ComBat, then recombine, then SVA (user design)

User: "since we know that sex and slide are confounded, couldn't we do ComBat
for slide and any other known technical variable for each sex strata, then
recombine the results to do SVA with PTSD, sex, age, cell proportions, and the
smoking proxy for any remaining unmeasured technical variance?"

Adopted. Rationale: cohort-wide, R2(sex ~ slide) = 0.874 -- five slides are
all-male, five all-female -- so a single ComBat pass must protect sex against a
variable nearly collinear with it. WITHIN a sex stratum that confounding
dissolves: slide becomes an ordinary batch factor. ComBat preserves the mean of
the matrix it is given, so correcting each stratum separately leaves
between-sex differences untouched by construction rather than by assumption.
PTSD remains balanced across slides inside both strata
(R2 = 0.057 female, 0.047 male), which is the condition that makes ComBat safe
for the exposure of interest.

### Within-stratum batch structure
  Female: n = 45, 7 slides, min batch 3, R2(PTSD ~ slide) = 0.057
  Male:   n = 42, 7 slides, min batch 1, R2(PTSD ~ slide) = 0.047
Only slides s6 (6F/1M) and s11 (3F/5M) carry both sexes.

### Pooling s6 + s11 was tested and REJECTED
The singleton male on s6 raised the question of merging s6 with s11. Testable,
because 9 females span both slides. Over the 30,000 most variable probes,
probe-centered slide centroid distances within females:
  s6 vs s11 = 152.2 -- rank 21 of 21 pairs, i.e. the MOST dissimilar pair
  s7 vs s8  =  72.9 (typical pair, permutation p = 0.59)
  median female pair ~ 90
s6 vs s11 separates on PC1 (R2 = 0.52, p = 0.029) and its centroid distance has
permutation p = 0.015 (2000 permutations). These two slides carry genuinely
different effects, so pooling them would hide a real batch difference inside one
term. Rejected on the evidence.

### mean.only = TRUE, set explicitly for BOTH strata
sva 3.46.0's ComBat silently forces mean.only = TRUE when any batch has a single
sample ("Note: one batch has only one sample, setting mean.only=TRUE"). Left
implicit, the male arm would receive a mean-only correction and the female arm a
mean-and-variance correction -- an undocumented asymmetry between strata. Set
explicitly for both so the two arms are treated identically. Slide differences
are thereby modeled as location shifts, which is the defensible component at
these batch sizes. No samples dropped; n stays 87 (45 F / 42 M).

### Array position: fixed covariate in the model, NOT a second ComBat pass
User challenge: "why is position being accounted for with SVA when we have that
measured? shouldn't that also be corrected for with ComBat since it's known?"

The first half of that is right and the original plan was wrong: leaving a
measured factor to SVA discards information. But ComBat is not the only way to
use a measured factor, and for position it is the wrong one.

Slide REQUIRES ComBat because it cannot enter the combined model at all: it is
87% collinear with sex (R2 = 0.874), so `~ ptsd + sex + slide` is effectively
unestimable. Correcting within strata is the only route.

Position has no such problem. 8 levels, 10-12 samples each, 4-6 per sex per
position, R2 vs sex = 0.014 and vs PTSD = 0.075. It is directly estimable as a
fixed effect in both the combined and the stratified arms, which propagates its
uncertainty into the standard errors instead of treating ComBat-adjusted values
as if they were observed data.

Three designs were run on the same matrix and scored on the 20,000 most variable
probes (lambda from the PTSD coefficient; residual R2 after fitting the design):

  arm                                   num.sv  lambda  slide  pos     df F/M
  A ComBat(slide), position -> SVA          15   1.071  0.0829  0.0786   20/17
  B ComBat(slide) + ComBat(position)        25   1.640  0.1058  0.0277   10/7
  C ComBat(slide), position as covariate    15   1.077  0.0914  0.0000   13/10

Arm B -- the "correct it since we measured it" reading -- fails three ways at
once. lambda inflates to 1.64, num.sv rises from 15 to 25, and residual slide R2
goes UP (0.0829 -> 0.1058): the second pass treats the already-corrected matrix
as raw data and partially undoes the first correction. It also produces 5 probes
at p < 1e-5 within a 20,000-probe slice where arms A and C produce none. This is
nygaard2016batch's anticonservatism measured directly rather than cited.

Arm C absorbs position completely (residual R2 = 0, by construction, since it is
in the design) and stays calibrated (lambda 1.077 vs 1.071). Its cost is 7
position dummies. ADOPTED: ComBat for slide within sex strata, array position as
a fixed covariate in the model, SVA for whatever remains unmeasured.

### Worker count: run with 8, document 4 (user ruling)
The tutorial's bash commands show `--workers 4`, which is a sensible default for
a reader's laptop. The recompute here is executed with `--workers 8` for speed.
Per user instruction the documented value stays at 4; only measured runtimes in
the comments are updated to match the actual run.

### Slide vs position vs well vs plate: which spatial factors are identifiable
User question: "There may only be 7 positions in this subset, but aren't there 96
positions possible on a full plate?"

Clarification first: there are 8 positions in the subset (R01C01..R08C01), all
present, 10-12 samples each. The "7" is the number of dummy COLUMNS 8 levels
contribute to the design matrix, not the number of levels.

The 96 refers to the plate WELL, and the distinction is the whole point. The
sample sheet is exactly one complete plate:

  12 slides x 8 positions = 96 wells = 96 samples, complete grid, 1 per well.

  factor    levels  samples/level  identifiable?
  slide         12        8         yes, but R2 = 0.874 with sex -> ComBat in strata
  position       8       12         YES - each chip row is reused on all 12 slides
  well          96        1         NO - 1:1 with sample; design would be saturated
  plate          1       96         NO - single plate, no between-plate contrast

array_pos is the ROW on the EPIC BeadChip (R01-R08), not the well. Rows are
reused across every slide, and that reuse is exactly what makes position
estimable as a fixed effect: 12 samples share R01C01, 12 share R02C01, etc.
Coding well identity instead would give 96 levels for 96 samples and zero
residual df.

Plate is NOT modeled because this dataset is a single plate. In a multi-plate
study plate would be a leading ComBat candidate, being the coarsest processing
batch. This is stated in the tutorial so students do not copy a single-plate
design into a multi-plate study unchanged.

### Are the SVs worth their degrees of freedom? (k = 0 test)

USER: "The SVs don't seem to explain much of the variance any more, are they still
worth including as covariates at all?"

A fair challenge, and the variance-explained table I had presented could not answer
it. Two reasons pve is the wrong metric here:

1. The 17.1% figure is a share of the **residual** variance — what remains after
   `ptsd + sex + age + smoke + position + 6 cell proportions`. 17% of the leftover
   is not a small structure to leave unmodeled.
2. A component can explain little variance and still bias a coefficient. What
   matters is whether it correlates with the exposure, not how much variance it
   carries.

So the design was tested directly with and without the SV block, full probe set
(756,251), adopted design otherwise identical:

| | no SVs | 6 SVs |
|---|---|---|
| Parameters | 18 | 24 |
| Residual df (combined) | 69 | 63 |
| **lambda** | **0.846** | **1.022** |
| Median SE of PTSD coefficient | 0.05577 | **0.04988** |
| Probes p < 1e-5 | 3 | 5 |
| FDR < 0.05 | 0 | 0 |
| Residual slide R2 | 0.0858 | 0.0828 |

cor(t) = 0.931; cor(-log10 p) = 0.857; top-20 overlap = 6/20; median |dt| = 0.232,
99.9th percentile 1.517.

**Findings.**

- **Dropping the SVs deflates the statistics to lambda = 0.846.** Unmodeled
  between-sample structure stays in the residual, inflates s^2, and shrinks every
  t-statistic toward zero. The chapter's own prose warns students that lambda << 1
  signals over-correction or underpowering; shipping a design that does this would
  contradict the lesson.
- **The SV block buys precision rather than costing it.** Six extra parameters at
  n = 87 should widen the SEs, but the median SE *falls* 11% and 77% of probes gain
  precision. The block absorbs ~18% of per-probe residual variance, and that
  reduction in s^2 outweighs the df loss. The usual "SVs cost df" tradeoff does not
  bind at this k.
- **Rankings are not robust to the choice.** Only 6 of the top 20 probes are shared
  between k = 0 and k = 6, so this decision changes which CpGs the tutorial presents.

**What the SVs are NOT doing:** residual slide R2 barely moves (0.0858 -> 0.0828).
ComBat and the position term already handle chip effects. The SVs remove unmeasured
structure with no name in this dataset — cell-composition estimation error, cryptic
relatedness, ancestry substructure, hybridization variability.

**Consequence for k.** k = 6 had been chosen on residual slide R2, which is flat
across k and therefore uninformative. Since additional SVs now look cheaper than
assumed, k is being re-selected on lambda and median SE across k = 6/8/10/15, taking
the first k SVs in SVA's native order (no reordering by chip loading).

USER ruling: "Re-test k=6/8/10/15 on lambda and SE and use SVA's native order"

**Superseded by this section:** the earlier argument that extra SVs "spend df on
variance that is neither slide, position, nor smoking." Unexplained residual
variance is worth removing whether or not it can be named.

### Choosing k on calibration and precision (final)

The user's ruling after the k = 0 test: re-test k = 6 / 8 / 10 / 15 on λ and median SE,
taking the first k SVs in SVA's native order (no reordering by chip loading).
`sweep_k_lambda.R` ran the full 756,251-probe limma at each k with the adopted design
(`~ ptsd + sex + age + smoke + pos + 6 cell proportions + k SVs`, stratified ComBat on
slide already applied). Results:

| k | params | resid df | λ | median SE | SE vs k=0 | p < 1e-5 | R²(smoke ~ SVs) | df F / M |
|---|---|---|---|---|---|---|---|---|
| 0 | 18 | 69 | 0.846 | 0.05577 | — | 3 | — | 28 / 25 |
| **6** | **24** | **63** | **1.022** | **0.04988** | **0.949** | **5** | **0.39** | **22 / 19** |
| 8 | 26 | 61 | 1.046 | 0.04938 | 0.936 | 7 | 0.44 | 20 / 17 |
| 10 | 28 | 59 | 1.072 | 0.04879 | 0.923 | 9 | 0.46 | 18 / 15 |
| 15 | 33 | 54 | 1.071 | 0.05025 | 0.945 | 7 | 0.64 | 13 / 10 |

**Adopted: k = 6.** Four reasons, in the order they bind:

1. **Calibration.** λ = 1.022 is the closest to 1 of any arm tested. k = 8 (1.046) is
   acceptable; k = 10 and k = 15 (1.07) drift into mild over-dispersion, and the p < 1e-5
   count rising monotonically with k (5 → 7 → 9) is the signature of variance being
   removed from the denominator faster than it is being removed from the signal.
2. **Precision costs almost nothing.** Median SE at k = 6 is 0.04988 against a minimum of
   0.04879 at k = 10 — k = 6 is 2.2% off the optimum while spending 4 fewer parameters.
   The SE curve is U-shaped: it bottoms at k = 10 and turns back up by k = 15, where more
   df are spent than variance is removed.
3. **The smoking proxy gets re-absorbed as k grows.** R² of the proxy on the SV block goes
   0.39 → 0.44 → 0.46 → 0.64. The proxy is in `mod0` precisely so SVA will not reconstruct
   it, but that protection is partial and it degrades with k. At k = 15 nearly two-thirds
   of the proxy lives inside the SV block, which double-counts smoking and makes the
   proxy's coefficient uninterpretable for students reading the model.
4. **The male stratum is the binding constraint.** The snakemake pipeline reuses the
   whole-cohort SVs inside each single-sex arm, so every SV costs a df there too. At
   k = 6 the male arm retains 19 residual df; at k = 15 only 10.

Residual slide R² is flat across k (0.083 / 0.083 / 0.084 / 0.091) — confirming again that
ComBat plus fixed-effect position, not the SVs, are what handle chip structure. Max R² of
any single SV with PTSD is stable at 0.040–0.041, so no arm is at risk of an SV proxying
the exposure. Per-stratum QR rank checks passed at every k.

Superseded by this section: the earlier k = 8 recommendation (chosen on the k-sweep table
before position entered the model as a fixed covariate) and both the ~80%-variance and
scree-plot selection rules, which do not apply to SVA's normalized, unordered components.

Figure: `data/05_k_selection.png`; table: `data/05_k_selection.csv`.
