# CHANGES — EWAS tutorial rebuild on Zhou EPIC `mask.cm` v8.1

Scope: chapters 00–08 of the Quarto tutorial in `repo/`, rebuilt so that the
executed code, the rendered numbers, and the prose all agree under the Zhou
**v8.1** recommended masking (read with **YAME v1.40**) instead of the legacy
`EPIC.hg38.population.manifest`. Dataset unchanged: Grady Trauma Project
(GSE132203), EPIC v1, n = 87 after QC (32 PTSD cases / 55 controls), ~95%
African American.

---

## 1. Masking: legacy manifest → v8.1 `mask.cm`

**What changed.** Probe filtering no longer reads
`EPIC.hg38.population.manifest`. It reads Zhou's v8.1 `EPIC.hg38.mask.cm`
(24 mask terms, 866,553 probes, ~360 kB packed) via the YAME binary, and
applies the **eight** codes appropriate to an admixed African-American cohort:

| Mask code | Probes flagged |
|---|---|
| `M_mapping` | 428 |
| `M_nonuniq` | 28,833 |
| `M_SNP_AFR_1pt` | 35,420 |
| `M_1baseSwitchSNP_AFR_1pt` | 445 |
| `M_2extBase_SNP_AFR_1pt` | 20,640 |
| `M_SNP_EUR_1pt` | 17,076 |
| `M_1baseSwitchSNP_EUR_1pt` | 259 |
| `M_2extBase_SNP_EUR_1pt` | 8,887 |
| **union of the eight** | **90,832** |

Cross-validated two ways: the YAME `unpack -a` route and a direct read of
`EPIC.hg38.mask.tsv.gz` both give union = 90,832, `identical() == TRUE`, with
all 90,832 inside the 865,918-probe tutorial universe.

**Numbers that moved.**

| Quantity | Legacy (as previously rendered) | v8.1 (now) |
|---|---|---|
| Probes masked | 124,132 (`MASK_general_AFR`) | 90,832 (8 codes) |
| Probes retained after filtering | 724,282 | **756,273** |
| Detection-p drops | 1,928 | 1,928 (unchanged) |
| Sex-chromosome drops | 19,627 | 19,627 (unchanged) |
| Retained if EUR flags only (counter-example) | 752,100 | 792,826 |

Net: **+31,991 probes (+4.42%)**. The direction is counter-intuitive and is now
called out explicitly in ch03: the legacy `MASK_general_<POP>` bundles repeat
and copy-number terms that v8.1 has no equivalent for, so the *older, single-
ancestry* mask is the more aggressive one.

**Legacy path retained in the text.** Per request, ch03 still documents that the
sub-population manifest exists (`EPIC.hg38.manifest.pop.tsv.gz`, 5 super-
populations / 26 sub-populations) and shows how to use it, with the
version-contrast table above so a reader can see the cost of the choice.

## 2. Full recompute of chapters 03–08

Every checkpoint downstream of the mask was regenerated rather than patched:
`03_grs_filtered.rds` (756,273 × 96) → `05_mvals_combat.rds` (756,251 × 87) →
`05_sva.rds` → `06_ewas.rds` → the Snakemake pipeline run → `08_annotation/`.

## 3. All probes tested, not a 20,000-probe subsample

**What changed.** The tutorial previously ran the EWAS on a seeded random slice
of 20,000 CpGs (`set.seed(42)`) for speed. That is gone. Both the in-notebook
EWAS (ch06) and the Snakemake pipeline (ch07) now run the **full matrix**.

**Numbers that moved.** Probes tested: 20,000 → **756,251**. Runtime claims were
re-measured, not estimated:

| Arm | Wall time |
|---|---|
| Combined | 1 h 54 m (113 m 58 s) |
| Female | 1 h 55 m (115 m 20 s) |
| Male | 1 h 03 m (62 m 59 s) |
| BACON | ~23 s per arm |
| METAL | ~90 s over 756,251 markers |

The "~15 seconds" / "~11 s" runtime claims and the `#| fig-cap` references to
"20,000 most variable probes" were replaced. A runtime callout with a 5-row cost
table was added to ch07, and the prose now says "roughly two hours on a busy
workstation."

## 4. The probe-count arithmetic, reconciled

Three different counts were rendering across chapters. Settled by direct
computation:

```
756,273   probes pass ch03 filtering
  −  20   smoking-proxy panel CpGs (excluded: circularity)
  −   2   probes with non-finite M-values (cg17759086, cg01801182;
          β pinned at 0 or 1 → log2(β/(1−β)) = ±∞)
= 756,251 probes tested
```

**Fixed:** the stale **756,271** in `05_batch_effects.qmd` (~line 298) and
`07_pipeline.qmd` (line 127) — both now read 756,273 → 756,251 with the two
deductions named. `03_probe_filtering.qmd` claimed "Every probe that survives
here is tested in the EWAS — there is no further subsetting downstream," which
was false; it now lists both deductions and states that the remaining sense in
which nothing is subsetted is that there is *no random subsampling*.

The `724,282` and `124,132` that still render in ch03 are **correct** — they are
the legacy row of the mask-generation contrast table, not stale counts.

## 5. Batch correction and covariate design (ch05)

**Settled design:** ComBat on `slide` **within each sex stratum**
(`mean.only = TRUE`, `par.prior = TRUE`), strata recombined, then array
**position as a fixed covariate** (not absorbed into SVA — it is measured), then
SVA with **k = 6**.

Rationale added to the text: `sex` and `slide` are confounded in this cohort, so
a single pooled ComBat on `slide` cannot be run without removing sex; stratifying
and recombining avoids that. Position is measured, so spending SVs on it is
wasteful.

**k = 6 chosen from a λ/SE sweep** (now a figure + table in ch05):

| k | λ | median SE |
|---|---|---|
| 0 | 0.846 | 0.0558 |
| **6** | **1.0217** | **0.0499** |
| 8 | 1.046 | — |
| 10 | 1.072 | — |
| 15 | 1.071 | — |

k = 6 recovers 76.5% of the attainable SE gain (se_ratio 0.949) at the λ closest
to 1. A scree-style truncation was considered and rejected — surrogate variables
are not ordered by variance explained the way PCs are, so a variance-cumulation
rule does not transfer; this is now explained in the text rather than left
implicit. The per-SV variance-explained figure (`05_sv_pve.png`) and the
"are SVs worth including at all" check (`05_sv_worth_it.png`) were added.

**PC1 of the batch PCA:** 17% → **15.4%** (now computed inline, not hardcoded).
Strongest associations: Neu R² = 0.79 with PC1; chip with PC3 (0.66) and PC5
(0.69); sex with PC5 (0.51); PTSD max R² = 0.03.

**Cell proportions:** the model previously used "5 of the 6" estimated
proportions; it now uses all **six** (`CD8T, CD4T, NK, Bcell, Mono, Neu`). Rank
was checked explicitly — design is full rank either way; VIF for `ptsd` is 1.224
in both, so the six-proportion version costs nothing in collinearity.

## 6. Smoking, as requested

**Added to preprocessing** (was absent). A smoking proxy is built as **PC1 of 20
replicated smoking-associated CpGs** (anchored on `cg05575921`, *AHRR*), oriented
so higher = heavier exposure; PC1 captures **53%** of panel variance. Cited:
`@joehanes2016smoking`, `@zeilinger2013smoking`, `@shenker2013epic`,
`@elliott2014smoking`, `@bollepalli2019epismoker`, `@bauer2015gpr15`.

**The 20 panel CpGs are excluded from the tested set** (`05_sva.rds$excluded_n`
= 20) — testing a probe against a score built from that probe is circular. This
is the source of the 20-probe deduction in §4.

`smoke` enters the final model as a covariate **and** the SVA null model, so SVA
does not re-discover smoking as a surrogate variable. Residual overlap is
reported honestly: R² between `smoke` and the SV basis is 0.394 at k = 6
(0.636 at k = 15).

**Final model formula** (ch06 and ch07 now agree):

```
~ ptsd + sex + age + smoke + pos + CD8T + CD4T + NK + Bcell + Mono + Neu + SV1..SV6
```

design_ncol = 24, residual df = 63, n = 87.

## 7. EWAS and BACON (ch06)

Re-run on the refiltered matrix with the model above.

| Quantity | Previously rendered | Now |
|---|---|---|
| Probes tested | 724,282 (then 20,000) | 756,251 |
| λ raw | 0.9908 | **1.0217** |
| λ after BACON | 1.0197 | **1.023279** |
| Top probe | cg25717994 / cg12509946 | **cg15434749** |
| Top BACON p | — | 2.96414e-08 |
| Bonferroni threshold | — | 6.611562e-08 |
| Bonferroni hits (BACON) | 0 | 1 |
| FDR < 0.05 (BACON) | 3 | 2 |
| p < 1e-5 | 5 | 22 |

`06_ewas.rds$tt` is 12 columns × 756,251 rows with 0 NA in `bacon.p`.

**RETRACTED.** An earlier draft added a `bacon.p.t` column and prose about a
"BACON *t*-reference." It was unnecessary and has been removed everywhere —
`bacon.p.t`, "t-referenced", "t-reference", and "normal reference" now return
zero hits across all rendered chapters.

## 8. Snakemake pipeline (ch07)

Pipeline inputs rebuilt for v8.1 and re-run (combined + female + male arms,
BACON per arm, METAL meta-analysis). Pipeline and in-notebook results agree:
`cor(estimate, logFC) = 1`, max absolute difference 8.88e-15.

**λ:** raw 1.021 / 0.967 / 1.070 (combined / F / M); after BACON 1.023 / 1.014 /
1.013. **Median SE:** 0.050 / 0.078 / 0.089.

**Cross-strata concordance — the honest reading, newly added.**
`cor(F_est, M_est) = −0.0043`; sign agreement **49.9%** (chance); heterogeneity
**5.4%** of probes below P = 0.05 against a 5% null; median I² = 0.

Nine CpGs clear Bonferroni in the meta-analysis, but **eight of the nine are
driven by one stratum alone**. That table is now in the chapter, together with an
explicit warning that "a small p-value from a two-study meta-analysis is not
evidence of replication." Only `cg21102739` is nominally significant in both arms
with matching sign. The closing caveats were rewritten to end on
"the honest summary of this run is a null result."

New figures: `#fig-meta-qq`, `#fig-qq-overlay`, `#fig-concordance`, plus
`07_pipeline_summary.csv` and the per-arm BACON diagnostic panels.

## 9. Functional annotation and enrichment (ch08)

Re-annotated against the Zhou v41/v36 gencode manifests on the new result set.

**Gene-set tests — nothing survives multiple testing:**

| Method | Sets tested | FDR < 0.05 | Best |
|---|---|---|---|
| `gometh` GO | 22,792 | **0** | GO:0070831 basement membrane assembly, P = 1.79e-05, FDR = 0.407 |
| `gometh` KEGG | 372 | **0** | hsa05217 basal cell carcinoma, P = 0.00101, FDR = 0.376 |
| `methylglm` | — | **0** | hsa04070 phosphatidylinositol signaling, p = 0.00472, padj = 0.448 |
| `methylRRA` (GSEA) | — | **0** | hsa04340 Hedgehog signaling, p = 0.000463, padj = 0.0976, NES = 1.20 |

**KYCG regulatory features:** 863 tested, 226 nominal p < 0.05 (vs ~43 expected),
**41 at FDR < 0.05**, all positive, log₂OR 0.361–1.282. Per knowledgebase:
TFBSconsensus 783/211/**35** (top `ING2`, log₂OR 0.55, FDR 0.00328);
HMconsensus 60/13/**5** (top `H3K36me2`, 0.95, FDR 0.0144);
chromHMM 15/1/**1** (`1_TssA`, 0.49, FDR 0.00328);
CGI 5/1/**0** (`Island`, 0.27, p = 0.00784, **FDR 0.0835**).

**Corrected claim.** The chapter previously closed by saying KYCG "recovers a
genuine **CpG-island enrichment** (FDR = 4 × 10⁻⁴)". It does not — the CGI
`Island` term does **not** clear FDR (0.0835). That takeaway was rewritten: all
three gene-set tools return nothing at FDR < 0.05 while KYCG returns 41 of 863
features, a difference driven as much by set size as by biology, whose honest
reading is a single broad enrichment in promoter-proximal, TF-bound chromatin
rather than 41 discoveries. Median significant-feature `nD` = 56,061 = **7.4%**
of the 756,251 universe — added so a reader can see the set-size effect.

**Knowledgebase versions now stated:** `KYCG.EPIC.CGI.20210713`,
`KYCG.EPIC.chromHMM.20211020`, `KYCG.EPIC.TFBSconsensus.20211013`,
`KYCG.EPIC.HMconsensus.20211013`. The earlier draft cited
`KYCG.EPIC.CGIposition.20211013`, which is not the database used.

KEGG pathway names were backfilled from `https://rest.kegg.jp/list/pathway/hsa`
(372 pathways) because methylGSA returns bare IDs.

**DMR analysis (comb-p).** At a defensible seed (1e-15) the run returns **0
regions**; the chapter shows a demonstration seed (1e-8) yielding 1 region and
labels it as a demonstration, not a finding.

A tutorial-caveat section was added before the closing checklist.

## 10. `data/CGI.png` and the region distances

The user-supplied figure was installed at `repo/data/CGI.png` (40,542 bytes;
referenced from `00_setup.qmd` line 27). No schematic was drawn.

Prose reconciled to the figure:

| Region | Was | Now |
|---|---|---|
| CpG island span | 500–1500 bp | 500–1500 bp (unchanged) |
| GC content | > 60% | > 60% (unchanged) |
| Shore | 0–2 kb | 0–2 kb (unchanged) |
| Shelf | 2 kb–4 kb | 2 kb–4 kb (unchanged) |
| Open sea | "greater than 2 kb" | **"greater than 4 kb"** |

Also corrected per user: the **70–80% methylated** figure refers to
**genome-wide CpGs**, not to CpG-island CpGs (islands are predominantly
*un*methylated). The sentence was rewritten accordingly.

## 11. Citations and residual consistency

- `joehanes2016smoking` was cited but undefined → added to `references.bib`
  (Circulation: Cardiovascular Genetics 9(5):436–447, 2016;
  doi 10.1161/CIRCGENETICS.116.001506).
- `aryee2014minfi` — now cited in `01_qc.qmd`.
- `smith2011grady` — now cited in `00b_dataset_catalog.qmd`.
- `hannum2013clock`, `horvath2013clock` — still defined but uncited (cosmetic;
  left in the bibliography deliberately).
- Undefined-citation sweep (`?@key`) across all rendered chapters: **zero**.
- `00b_dataset_catalog.qmd`: the shortlist-funnel figure
  (`data/00b_shortlist_funnel.png`) was referenced but absent → the figure block
  was removed and the shortlist described in prose (70 studies screened → 20
  shortlisted → GSE132203).
- `05_batch_effects.qmd` had a duplicated `n.sv` display → fixed.
- ch08 lines 79 / 107 / 109 / 110 use `CGIposition` as a **v36 manifest column
  name** — legitimate, deliberately left alone.

## 12. Verification harness

`verify_txt/` is generated from `repo/_site/*.html` by stripping base64 payloads
(`data:…;base64,…` → `[B64]`), `<script>` → `[JS]`, and `<style>` → `[CSS]`.
Stale-token audits **must** run against `verify_txt/`, not raw HTML — embedded
figure payloads produce false positives (e.g. the strings `20k`, `4.5%`, `0.77`
appear inside base64 blobs while the source `.qmd` contains none of them).

Tokens swept to zero: `bacon\.p\.t`, `t-referenced`, `t-reference`,
`normal reference`, `724,282`/`724282` (outside the legacy-contrast table),
`cg12509946`, `CGIposition\.20211013`, `HNF1A`, `0, 1, and 7`, `4\.1 × 10`,
`4 × 10⁻⁴`, `756,271`, `756,273` (as a tested count), `random 20k`,
`random slice`.

Tokens verified present: `756,251`, `1.0217`/`1.022`, `1.023`, `cg15434749`,
`49.9`, `5.4%`, `−0.004`, `41 of 863`, `CGI.20210713`, `v8.1`, `1 h 54 m`,
`2.0 × 10`, `0.084`.

## 13. Known-benign render noise

Left unfixed deliberately: `fonts.googleapis.com` and `cdnjs.cloudflare.com`
(MathJax polyfill) return 403 through the sandbox proxy; the pages render
correctly without them. Also benign: the bubblewrap `--bind-fd` warning and
`Fontconfig warning: using without calling FcInit()`.

## 14. Not re-run

- **Chapter 04** (cell composition) was not re-executed: `FlowSorted.Blood.EPIC`
  is not installed in the `methyl` environment. The existing cell-proportion
  estimates (`04_cc_full.rds`) are reused, and they are upstream of the mask
  change, so the v8.1 refilter does not invalidate them.
- `recompute_06_final.R` has not been re-executed since its `saveRDS` patch;
  `06_ewas.rds` currently carries columns written by a manual `match()` backfill.
  The values were verified against the pipeline output (§8) and agree to 8.9e-15.

---

## 15. Sample scope: where 96 becomes 87

**The question.** Nine of the 96 QC-passing arrays have no recorded PTSD status.
Earlier drafts collapsed this into a single "n = 87 after QC" sentence, which
misattributes a metadata gap to a quality failure and hides the fact that the
count changes exactly once, at the step where the exposure enters.

**Chapter 01** now separates *missing values in variables we intend to model*
from *variables absent for every sample*. Two new chunks compute the accounting
live rather than asserting it:

| variable | missing of 96 |
|---|---|
| `ptsd` | **9** |
| `sex` | 0 |
| `age` | 0 |
| `childhood_abuse` | 2 |
| `mergedcapsandpsswinthin30days` | **9** |

The CAPS/PSS severity score is missing for the *same* nine samples, so no
fallback outcome exists. The nine are 4 female / 5 male, ages 22–65, spread
across 8 of 12 chips — incidental, not structured, which is the missing-at-random
check the chapter now walks through. A `callout-important` states the
distinction: samples dropped for *quality* are samples you do not trust, while
samples dropped for *missing phenotype* are samples you trust but cannot model.

**Which steps run on which count**, now stated explicitly in ch04, ch05 and ch06:

| step | n | why |
|---|---|---|
| QC (ch01) | 96 | all pass |
| Normalization (ch02) | 96 | no model |
| Probe filtering (ch03) | 96 | no model |
| Cell deconvolution (ch04) | 96 | `estimateCellCounts2` takes no model |
| PCA diagnostics (ch05 §2) | 96 | descriptive; no model |
| Smoking proxy (ch05 §5) | 87 | captures variance of the set the EWAS fits |
| Stratified ComBat (ch05 §4) | 87 | takes a protection model containing PTSD |
| SVA (ch05 §6) | 87 | `sva()` requires the full model; see below |
| EWAS (ch06, ch07) | 87 | 32 cases / 55 controls |

**Why SVA cannot be moved to 96.** `sva()` requires both a full model and a null
model, and the full model is what prevents the surrogate variables from absorbing
the exposure. Because `ptsd` is `NA` for nine samples, `model.matrix` returns 87
rows — it drops incomplete rows silently — so a call against a 96-column matrix
fails with `non-conformable arrays`. Forcing 96 through would require either
removing `ptsd` from the full model, which makes it identical to the null model
so that nothing is protected and the SVs are free to regress away the signal
under test, or imputing the exposure. The 87 gate is where SVA's definition puts
it, not a preference.

**Consequence worth noting** (now in ch05 §4). Chip `201114400024` carries 6
females and **2** males across all 96 arrays, but 6 females and **1** male among
the 87. One of the two males is among the nine. That singleton is what forces
`mean.only = TRUE` in the male stratum — so the incomplete phenotype table did
not just cost nine samples, it changed which batch corrections are estimable.
ComBat *could* have been run on 96, where the minimum male chip is 2 and variance
scaling is identifiable; this was considered and declined, because ComBat's
protection model also contains PTSD, the matrix still has to narrow to 87 before
SVA regardless, and the singleton is more valuable as a worked example of missing
metadata propagating into a methods choice.

## 16. Prose revision pass

The author revised prose across all eleven source files and the README. Those
revisions were installed verbatim; previous copies are preserved in
`repo_backup_pre_userprose/` and `repo_backup_pre_userprose2/`. Substantive
changes made on top of the author's text, all of them corrections:

- **ch05** — `mean.only = TRUE` was described as correcting "each chip's mean and
  variance". Corrected to "mean but not its variance", which is what the argument
  does and the reason the male/female asymmetry matters.
- **ch05** — the "Where 96 becomes 87" callout linked `04_cell_composition.qmd`
  for a sample-funnel section that had been deleted. Repointed to
  `01_qc.qmd#what-the-phenotype-file-is-missing`.
- **ch05** — the callout named only SVA as running on 87, leaving the impression
  ComBat ran on 96. It does not: `05_mvals_combat.rds` is 756,251 × 87.

Deleted by the author and **not reinstated**: the ch04 sample-funnel section, the
ch05 `chip-position-grid` chunk, and the "one full plate" paragraph.

## 17. US spelling throughout

All source files, both `repo/` and the packaged `gh/EWAS-tutorial/` copies, plus
all 14 recompute and figure scripts, were converted to US spelling — 42 files.
Word-boundary matched, then re-scanned to zero residuals.

| category | tokens |
|---|---|
| `-ise` → `-ize` | analyse/analysed/analysing, meta-analysed, summarise, residualise, randomisation |
| `-ll-` → `-l-` | modelled, modelling, unmodelled, labelled, signalling |
| `-our` → `-or` | colour (32×), coloured, behaviour, neighbour(s)(ing) |
| `-re` → `-er` | centre, centred |
| `-gue`/`-ce` | catalogue (21×), catalogued, analogue, licence |
| misc | artefact(s), ageing, grey → gray |

Code identifiers were converted deliberately: ggplot2 accepts `color`/`colour` as
synonyms and R accepts both `grey`/`gray` colour names, so `colour =` → `color =`,
`scale_colour_manual` → `scale_color_manual`, and `"grey45"` → `"gray45"` are
behaviour-preserving. All 25 R scripts were re-parsed after the substitution to
confirm no syntax errors, and the regenerated figures are byte-comparable in
appearance.

Left unconverted, correctly: `analyses` as the plural of *analysis* (US-correct);
`meta_analysis` / `dmr_analysis` config keys and the
`PTSD_ewas_meta_analysis_results_1.txt` filename (already US spelling); and
`characteristics` where it names GEO's own series-matrix field.

## 18. Final render, verification pass, and publishing

### Full-site re-render

The whole book was re-rendered from a clean state after the US-spelling pass
(`bash qrender.sh render`), producing eleven HTML files in `_site/`. Two earlier
launch attempts had died silently; the cause was environmental rather than
analytical, and is recorded in the tooling notes rather than here.

### Rendered-prose audit

Rendered text was extracted from every chapter (tags and base64 payloads
stripped) and scanned rather than trusting the sources. Results:

| check | outcome |
| --- | --- |
| British-spelling residuals | none (only `analyses`, the US-correct plural) |
| `756,251` (CpGs tested) | present in 23 places |
| `724,282` (superseded retained count) | absent from prose; appears only as the 4 legitimate legacy-mask contrast values |
| `1.0104` / `1.0198` (superseded λ) | absent |
| `715,992` (rejected union variant) | absent |
| `n = 87` / `n = 96` | both present, each in its correct scope |

The λ values that survive in the rendered text are the current ones
(λ_raw ≈ 1.02), rendered from the checkpoint via inline `sprintf` rather than
hard-coded, so they cannot drift from the data again.

### Two stale Git-LFS pointer files removed

`data/06_ewas_bacon_toptable.csv` and `data/06_ewas_toptable.csv` were not data
files at all — each was a 133–134 byte Git-LFS pointer stub left over from the
original upstream repository, naming an object that was never fetched into this
workspace. Nothing in the tutorial referenced them (every chapter reads the
`.csv.gz` form or the `.rds`), and `06_ewas_bacon_toptable.csv.gz` (60.8 MB,
regenerated under the v8.1 mask) is the real artifact. Both stubs were deleted
from `repo/data/` and from the packaged tree so that a reader cloning the
repository without `git-lfs` does not find a file that looks like results and
is not.

### Unverifiable comment softened

The `listDBGroups("EPIC")` code comment claimed "19 EPIC knowledgebases". The
count could not be re-verified — `knowYourCG` is not installed in the render
environment, and the chapter's enrichment chunks are `eval: false` for exactly
that reason — so the specific number was removed and the comment now lists the
knowledgebase categories without asserting a count. The verified numbers from
the recorded run (863 features tested, 41 at FDR < 0.05, across four
knowledgebases) are unchanged, because those came from the run itself.

### Citation audit

All five previously-flagged references (`hannum2013clock`, `horvath2013clock`,
`aryee2014minfi`, `smith2011grady`, `joehanes2016smoking`) are now both defined
in `references.bib` and cited in the text. No undefined citation keys remain.

### Repository name

The GitHub repository is `krferrier/Methylation-EWAS-tutorial` (hyphenated).
Thirteen references across `README.md`, `CITATION.cff`, `PUBLISHING.md` and
`zenodo.json` were updated: the Pages URL, the clone URL, the `cd` lines, the
`gh repo create` invocation, the CFF `repository-code` and `url` fields, and the
Zenodo description and related-identifier URLs. The Pages URL is
`https://krferrier.github.io/Methylation-EWAS-tutorial/`. Note that the on-disk
staging directory retains its original name; only the references were renamed.

## 19. Site theme, figure restyle, and the ch08 DMR reconciliation

This section covers the final presentation pass and the last set of numbers that
moved. It is the only pass in which published *figures* changed; the underlying
analysis was already frozen and no model was re-fit.

### A compiled theme replaces the bare Bootstrap preset

`tutorial/_quarto.yml` previously set `theme: cosmo`. It now sets a two-element
theme stack, `[cosmo, theme.scss]`, so the site is Cosmo with a compiled SCSS
layer on top, and adds `fig-dpi: 200` so raster figures embed at the same
resolution the restyle script writes them at.

`tutorial/theme.scss` is new (326 lines, 8,545 bytes). It defines a *Scholarly
teal + warm gray* palette, sets `$font-family-sans-serif` to Source Sans 3,
`$font-family-serif` to Source Serif 4 (also used for headings via
`$headings-font-family`), and `$font-family-monospace` to JetBrains Mono, each
with a full platform fallback chain. A `@font-face` loop at line 100 emits
`src: url("assets/fonts/#{$file}.woff2") format("woff2")` for every face.

### Fonts are self-hosted, not fetched from Google

`tutorial/assets/fonts/` is new and carries nine woff2 faces (Source Sans 3
Regular/Italic/Semibold/Bold, Source Serif 4 Regular/Semibold/Bold, JetBrains
Mono Regular/Bold). Quarto copies them into the rendered site, and all eleven
pages in `docs/` were verified to carry two `@font-face` rules and to reference
only local asset paths — no page loads a webfont from a third party.

`tutorial/assets/fonts/otf/` additionally ships four OTF sources
(SourceSans3 Regular/Italic/Semibold/SemiboldItalic, 1.15 MB). These are *not*
redundant with the woff2: the browser consumes woff2, while `_setup.R` registers
the OTF files with `systemfonts` (lines 43-61) so that R figure text in a fresh
clone renders in the same typeface as the surrounding page. If registration
fails — no `systemfonts`, or a version too old — `.ewas_family` falls back to
`""` and plots use the device default rather than erroring. The whole assets tree
is 2.2 MB and is tracked; `.gitignore` deliberately does not exclude it, because
a clone without the fonts renders in the wrong typeface with no warning.

### Palette literals removed from chunk bodies

`_setup.R` now exports a named palette and a theme:

- `ewas_col` — `teal_dark #0F3D43`, `teal #1A6B75`, `teal_light #2A8F9B`,
  `sand #B8873F`, `plum #8C3A4A`, `gray_dark #3A362F`, `gray #6E675B`,
  `gray_light #B5AFA4`
- `pheno_pal` — `Control`/`Male`/`neg` on teal, `Case`/`pos` on plum,
  `Female` on sand
- `theme_ewas()`, installed globally at line 96 via `theme_set(theme_ewas())`

Thirteen hard-coded seaborn hexes were replaced with palette lookups across three
chapters: **eight** in `01_qc.qmd` (detection-p scatter, intensity scatter, sex
scatter, SNP-fingerprint colour ramp), **three** in `02_normalization.qmd`
(Type I/II density facets, β→M transform curve), **two** in
`03_probe_filtering.qmd` (ancestry-mask contrast bars). `"#4C72B0"` became
`ewas_col[["teal"]]`, `"#C44E52"` became `ewas_col[["plum"]]`, `"#DD8452"`
became `ewas_col[["sand"]]`. The only hex literal remaining in any `.qmd` is
`00b_dataset_catalog.qmd:80`, `styleEqual("GSE132203", "#FFF3CD")` — an
intentional highlight in a DT table, not a plot colour.

All fifteen referenced PNGs were regenerated through a single script,
`scripts/restyle_figs.R` (541 lines, 30,659 bytes), rather than by re-running
chapter chunks; that keeps the figure pass independent of the analysis
checkpoints. Every regenerated figure was inspected before acceptance. Two
figures needed layout work beyond the palette swap: the SVA bar panel and the
DMR locus panel, the latter accepted only at the fourth iteration with
`plot.margin = margin(5.5, 22, 5.5, 5.5)` and `breaks_width(50)` on the
coordinate axis. The four `bacon` diagnostic panels are JPGs written by
`bacon`'s own plotting methods and were deliberately left un-restyled.

### `index.qmd`: two removals

The dataset-rationale table lost its "Deposited epigenetic-age acceleration"
row, and the "One dataset, on purpose" callout was removed entirely. The
tutorial never validates an epigenetic clock, so advertising the answer key for
one promised a section that does not exist.

### The ch08 DMR numbers were wrong and are now reconciled

This is the substantive correction in this pass. The chapter had been describing
a stale comb-p run (`data/08_annotation/dmr_demo/`) while the v8.1 rebuild wrote
its output to `data/08_annotation/dmr/`. The stale directory has been deleted,
`run_combp.sh` was patched to point at the authoritative path, and the demo
command in the chapter now reads `-p data/08_annotation/dmr/PTSD_dmr_demo`.

Numbers that moved, all from the stale run to the v8.1 run:

| quantity | was | now |
|---|---|---|
| Šidák region p | 1.8 × 10⁻⁹ | **1.2 × 10⁻¹⁰** |
| best raw BACON p among the six CpGs | 1.5 × 10⁻⁴ | **1.7 × 10⁻⁵** |
| per-CpG −log₁₀ p range | 3.3 – 3.8 | **2.8 – 4.8** |

The Šidák value appears in three places — the region table, the `fig-dmr`
caption, and the chapter summary — and all three were updated. The region itself
is unchanged: `chr6:28,633,534–28,633,601`, 67 bp, 6 CpGs, all
`Promoter`/`N_Shore`, `genesUniq = ENSG00000271440;ENSG00000287279`, all
hypomethylated.

Two column headers were also wrong. The region table's second data column was
labelled "min single-CpG p" but holds the best raw BACON p *within the region*,
not the minimum over all CpGs; it is now labelled accordingly. And the
`dmr_cpg_display.csv` table caption now states the Δβ range explicitly
(**−3.5 to −5.6 percentage points**, all in the same direction) instead of
gesturing at "negative effect on the β-scale". Both CSVs behind those tables
were regenerated: `dmr_cpg_display.csv` (6 rows, 409 bytes) and
`dmr_cpg_annotation.csv` (7 rows, 650 bytes).

### `min_p` semantics corrected

The chapter previously said comb-p's `min_p` is "the Stouffer–Liptak–Kechris–
smoothed p-value". That is half right and misleading in the half that matters:
`min_p` is the smallest SLK-smoothed **FDR-corrected regional q** among the
region's CpGs. The raw SLK p at those same CpGs goes down to 6.0 × 10⁻¹⁹, so the
text now distinguishes the two — regional p 6.0 × 10⁻¹⁹ *before* FDR correction,
q 2.3 × 10⁻¹³ *after* — in both the sentence about the strongest local cluster
and the dedicated `min_p` subsection. Without that distinction a reader
comparing 2.3 × 10⁻¹³ against a probe-level p is comparing two different
statistics on two different scales.

### Column-name correction: `bacon.p.t` does not exist

An earlier decision recorded `bacon.p.t` as the p-value column in
`annotated.rds`. It is not a column in that object. The correct name is
**`bacon.p`**, and the DMR figure code and the comb-p invocation both use
`pvals = bacon.p`. This supersedes the earlier note.

### Two captions corrected against their own figures

- **SVA bars.** The caption claimed "the strongest chip-linked SV reaches
  R² ≈ 0.9". Under the k = 6 checkpoint the strongest is **SV5 at R² = 0.64**,
  and sex is entangled with chip in the same SVs, so the caption now names all
  three tracks (chip, sex, PTSD) and gives the largest PTSD association,
  **R² = 0.04**, as evidence the SVs stay near-orthogonal to the exposure.
- **KYCG enrichment.** The caption said the twelve features came "across four
  EPIC knowledgebases". Four were *tested*; only three are represented among the
  significant features, and the CGI knowledgebase contributes none. The caption
  now says exactly that, which also makes the "no CpG-island-relation feature
  appears" observation that follows it coherent rather than contradictory.

### Data-tier and record metadata

`data/08_annotation/dmr_demo/` was removed from the H_annotation tier and the
regenerated comb-p bundle replaced it, so the tier grew from 305 MB to
**306 MB** and its checksum changed to
`d23a810f60b791e582fd8a952ec907fcd9e8dedb37847b3abe1d689a2d02a8cf` in both
`MANIFEST.md` and `SHA256SUMS.txt`. The stale
`PTSD_ewas_annotated_results.csv.gz` was dropped from the repository data
directory (the `.bed` and the `_zhou.csv.gz` remain, and the `.csv.gz` was
redundant with them).

`get_data.sh` no longer ships a `REPLACE_WITH_RECORD_ID` placeholder and the
guard clause that aborted on it is gone; `ZENODO_RECORD` now defaults to
**22135216** and is still overridable from the environment for pinning an older
version. `MANIFEST.md` gained a header paragraph naming the concept DOI
**10.5281/zenodo.22135215** as the citable identifier, with the version DOI
**10.5281/zenodo.22135216** for this build. `zenodo.json` gained the creator's
affiliation (University of Colorado Anschutz Medical Campus) and ORCID
(0000-0002-8813-6871), and the file now ends with a newline.

### Known record-side defect (not fixable from any repository file)

The **published** Zenodo record's `related_identifiers` carry
`isSupplementTo: https://github.com/krferrier/EWAS-tutorial`, which does not
resolve — the repository was renamed to `Methylation-EWAS-tutorial` after that
metadata was drafted and Zenodo does not re-read `zenodo.json` after
publication. `zenodo.json` in this repository already has the correct URL, so
any future version will be right. Fixing the existing record requires editing it
in the web interface: open `https://zenodo.org/records/22135216`, click
**Edit**, change the `isSupplementTo` URL to
`https://github.com/krferrier/Methylation-EWAS-tutorial`, and **Publish**. This
is a metadata-only edit; it does not mint a new version and does not change
either DOI. Documented in `PUBLISHING.md` §"Known record defect".

### Packaging

`sync_gh_repo.sh` grew from 912 to **1,237 bytes** with a new step 1b that
copies `repo/theme.scss` into the packaged tree and `rsync -a --delete`s
`repo/assets/` across, so the theme and fonts travel with the sources instead of
being reconstructed by hand at publish time.

### Render verification

Full-site render exited 0 and wrote eleven HTML files. Every page was checked
for three things: at least one theme hex (`0F3D43` or `1A6B75`) present, two
`@font-face` rules present, and **zero** occurrences of the retired seaborn
hexes (`4C72B0`, `C44E52`, `DD8452`). All eleven pass all three. The rendered
`08_annotation.html` was additionally checked to contain each corrected number
(1.2 × 10⁻¹⁰, 1.7 × 10⁻⁵, 2.8–4.8, −3.5 to −5.6, 6.0 × 10⁻¹⁹, 2.3 × 10⁻¹³) and
*not* to contain the superseded ones (1.8 × 10⁻⁹, 1.5 × 10⁻⁴, `dmr_demo`,
`bacon.p.t`).

Two blocked outbound requests appear in the render log and are benign:
`fonts.googleapis.com` and `cdnjs.cloudflare.com/polyfill/v3/polyfill.min.js`,
both 403 via the sandbox proxy. Neither is referenced by any source file in this
repository — both come from Quarto's stock Bootstrap template — and with the
fonts now self-hosted the Google request has nothing to fetch that the site
needs.

---

## 20. Chapter 07 prose reconciliation, reference normalization, and the ch06 checkpoint restoration

### Chapter 07 prose reconciled against the author's revision

`07_pipeline.qmd` received nine prose edits supplied by the author, verified
line-for-line against the submitted source with a unified diff. The substantive
changes soften three claims that the data do not carry. The limma-vs-`glm`
comparison no longer asserts that "moderation buys stability at small *n*" and
instead notes that the ~22% standard-error discrepancy makes it important to be
transparent about which method produced a given estimate. The agreement between
the two engines is described as giving confidence in consistency rather than as
"exactly what a successful port looks like". The stratified-λ discussion now
says the smaller stratum is the more likely to show in- or deflation, rather
than instructing the reader which λ to "look at hardest". The sex-discordance
passage no longer calls the pattern "the single most important thing to notice
on this page"; it states that a small p-value with a discordant direction may
indicate a true sex difference in effect, and a new closing paragraph warns that
the per-stratum sample sizes are small for an EWAS and that any sex-specific
signal needs replication in a larger subset or the full cohort before it is
believed.

### Three typographical corrections

In `07_pipeline.qmd`: `on an remote` → `on a remote`; `The sample sized of each
strata` → `The sample sizes of each strata`; `statstically` → `statistically`.

### Cross-reference style normalized to "chapter NN"

The tutorial referred to its own chapters three different ways — `P4`,
`notebook 06`, and `chapter 4`. All are now `chapter NN`, zero-padded, matching
the file-name convention. Edits landed in `06_ewas.qmd` (nine references),
`07_pipeline.qmd` (two), and `08_annotation.qmd` (two). A tree-wide grep for
`\bP[0-9]+\b` and `notebook [0-9]+` across all eleven `.qmd` files now returns
nothing. Bracketed prose links of the form `[Setup](00_setup.qmd)` were left
alone: those are titles, not numbered cross-references.

### Two render-tree checkpoints were missing and one was stale

The v8.1 refilter had been carried out, but two intermediate checkpoints that
chapters 05 and 06 load at render time had never been regenerated under the new
mask, and `06_ewas.rds` on disk still held the pre-v8.1 fit. This did not show
up in the published HTML because `repo/data/06_bacon_summary.rds` — which is
what most of chapter 06's inline expressions actually read — had been written by
`scripts/recompute_06_final.R` under the correct mask. The published site was
therefore right; the reproduction tree was not.

Two scripts restore it, both taking every parameter from the frozen checkpoints
rather than re-deriving it:

- **`scripts/rebuild_checkpoints_v81.R`** (3,201 bytes) writes
  `data/03_grs_filtered.rds` as `02_funnorm_grs.rds[funl$retained_probes, ]` →
  **756,273 × 96**, which matches `03_filter_funnel.rds$retained` exactly; and
  `data/05_mvals_combat.rds` as stratified ComBat on slide (`mean.only = TRUE`,
  `par.prior = TRUE`, run separately within each sex, slide levels renamed
  `s1..sN`) with `cb_mod = model.matrix(~ ptsd + age + smoke + Neu + NK + CD4T +
  CD8T + Bcell + Mono)`. The `keep`, `mdk`, `propk`, `smoke`, and `k_selected`
  fields come from the existing `05_sva.rds`, which was opened read-only and not
  rewritten. The 20 smoking-panel CpGs in `05_smoking_proxy.rds$present` are
  excluded, leaving **756,251** probes tested — equal to `sva_o$tested_probes`.
- **`scripts/rebuild_06_rds.R`** (4,752 bytes) is the EWAS-and-BACON half of
  `recompute_06_final.R` with two side effects deliberately removed: it does not
  rewrite `06_bacon_summary.rds` (the on-disk copy carries extra fields chapter
  06 reads) and it does not redraw the Manhattan, volcano, or BACON diagnostic
  images (the on-disk copies are the restyled published ones). Design is
  `model.matrix(~ ptsd + sex + age + smoke + pos + CD8T + CD4T + NK + Bcell +
  Mono + Neu + SV)` with `ptsd` releveled to reference `"Control"`, coefficient
  `ptsdCase`, `SV = sva$SV[, 1:6]`, and `pos = factor(mdk$array_pos)`. A
  verification loop compares sixteen fields of the rebuilt object against
  `06_bacon_summary.rds`; all sixteen match.

The restored numbers, all of which agree with the previously published HTML:

| Quantity | Value |
|---|---|
| Design columns | 24 |
| Residual df | 63 (= 87 − 24) |
| Probes tested | 756,251 |
| λ before BACON | 1.022 |
| λ after BACON | 1.023 |
| Top CpG by BACON *p* | cg15434749 |
| Bonferroni-significant | 1 |
| FDR < 0.05 | 2 |

The design reconciles term by term: 1 intercept + ptsd + sex + age + smoke +
7 position + 6 cell + 6 SV = 24.

Pre-restoration copies of `03_grs_filtered.rds` and `06_ewas.rds` were kept
outside the repository tree during the rebuild and are not shipped.

### Full-site re-render

All eleven chapters were re-rendered from the restored checkpoints. The render
log shows `[ 1/11]` through `[11/11]` and no errors; Quarto emits one
`Output created` line for a website project, not one per page. Eleven HTML files
were written.

The rendered `06_ewas.html` was checked to contain each anchor above and *not*
to contain the superseded values (`756,271`, λ `1.010`, a 25-column design). The
CpG `cg25717994`, which an earlier stale fit had put at rank 1, now appears at
rank 3 of the top-six table with BACON *p* = 2.83 × 10⁻⁷ — the position the
published numbers give it. A tree-wide scan for the pre-v8.1 retained-probe
count `724,282` returns zero occurrences in all eleven pages.

The two blocked outbound requests documented in §19 (`fonts.googleapis.com` and
the `cdnjs.cloudflare.com` polyfill, both 403 via the sandbox proxy) appear
again and remain benign for the same reasons.

## 21. The comb-p seed threshold, and two upstream staleness bugs it uncovered

You reported that the DMR seed p-value in the pipeline's `config.yml` was
`1e-15` rather than the intended default `1e-4`, and asked for the comb-p run
and the surrounding prose to be redone. Fixing the seed alone would have
published a spurious result, because two separate stale inputs were sitting
underneath the DMR section. All three are fixed here.

### 21.1 The seed threshold

`--seed` is the p-value a CpG must reach before comb-p will *start* a region
there. It is not a significance threshold for the output — the Šidák-corrected
region p is. Setting it at genome-wide-significance levels prevents any region
from ever being seeded.

| seed | regions found | regions passing Šidák < 0.05 |
|---|---|---|
| 1 × 10⁻¹⁵ | 0 | 0 |
| 1 × 10⁻⁴ (default, now used) | 3 | 1 |

The strongest single CpG in this dataset reaches p = 3.0 × 10⁻⁸, so nothing
could ever seed at 1 × 10⁻¹⁵. A new `.callout-important` in chapter 08 ("Set the
seed threshold on purpose") documents the distinction and the empty-result
symptom.

### 21.2 Stale BED input

The `make-bed` chunk was reading `PTSD_ewas_annotated_zhou.csv.gz` — the July
pre-v8.1 table, 724,282 probes with pre-v8.1 p-values — rather than the v8.1
refit in `06_ewas.rds` (756,251 probes). Run at the corrected seed on the stale
BED, comb-p returns **two** regions passing Šidák < 0.05: the chr6 MHC region
plus `chr16:3112431-3112453` at Šidák 0.019. On the v8.1 BED that second region
does not exist. It was an artifact of the superseded p-values.

The chunk now reads `readRDS("data/06_ewas.rds")$tt` joined to the YAME v8.1
`EPIC.ordering.tsv.gz` + `EPIC.hg38.coord.tsv.gz`, filtered to primary contigs.
That writes 756,234 of the 756,251 tested CpGs; the 17 dropped fall on unplaced
or alt contigs. This was cross-validated against an independent coordinate route
(the GENCODE v41 manifest, which also carries hg38 coordinates keyed on
`probeID`): the v41 route yields 756,251 rows, and the 17-row difference is
exactly the non-primary contigs, enumerated one by one.

### 21.3 Stale annotation layer (`annotated.rds`)

A third staleness, found while auditing chapter 08's data references. The
shipped `data/08_annotation/annotated.rds` was a 724,282-row table on the **hg19
Illumina** schema (`UCSC_RefGene_Name`, `Relation_to_Island`, `bacon.pval`),
dated Jul 30 — it predated the entire Zhou hg38 rebuild. The v8.1 replacement
(756,251 × 22) had been generated on Aug 26 and was then overwritten by a
restore from the July render bundle in `~/Downloads`.

This did **not** affect any published number: every chapter-08 chunk that reads
`annotated.rds` or `PTSD_ewas_annotated_zhou.csv.gz` is `eval: false`, and the
displayed tables and figures come from precomputed CSVs and PNGs whose values
were confirmed to match `06_ewas.rds` exactly. The stale file was a shipped
data artifact that disagreed with the prose describing it, not an input to the
render.

The annotation layer was rebuilt through the code path chapter 08 documents —
`fread("data/06_ewas_bacon_toptable.csv.gz")`, joined against the GENCODE v41
and v36 manifests. Every published number reproduced:

| quantity | published | rebuilt |
|---|---|---|
| annotated CpGs | 756,251 | 756,251 |
| gene-mapped | 647,034 (85.6%) | 647,034 (85.6%) |
| BIOS eQTM-annotated | 10,803 (1.4%) | 10,803 (1.4%) |
| Bonferroni threshold | 6.61 × 10⁻⁸ | 6.612 × 10⁻⁸ |
| min BACON *p* | 2.96 × 10⁻⁸ | 2.964 × 10⁻⁸ |
| Bonferroni / FDR < 0.05 | 1 / 2 | 1 / 2 |
| promoters in top ten | 4 | 4 |
| rank of cg16340178 | 8 | 8 |

`top10_display.csv` and `eqtm_display.csv` came back **byte-identical** to the
shipped files. `eqtm_top_hits.csv` differs only in the last one or two decimal
digits of some floats (max absolute difference 1 × 10⁻¹⁵) — CSV round-trip
precision, with the probe column identical. `08_feature_distribution.png` and
`08_island_distribution.png` were left untouched; their percentages were
confirmed correct against the rebuild.

The two GENCODE manifests needed for this had been reported unavailable earlier
in the project. They fetch cleanly from
`github.com/zhou-lab/InfiniumAnnotationData` (v41: 62,477,035 B, 866,554 rows;
v36: 29,916,134 B, 865,919 rows).

### 21.4 The DMR that survives

Rerun at seed 1 × 10⁻⁴ on the v8.1 BED, comb-p reports three candidate regions
and one survives region filtering:

| region | CpGs | Šidák *p* |
|---|---|---|
| chr6:28,633,493–28,633,701 | 11 | **2.6 × 10⁻¹⁶** |
| chr3:138,608,544–138,608,573 | 2 | 0.45 |
| chr20:32,018,064–32,018,066 | 1 | 1.00 |

Autocorrelation in the input: correlation 0.054 out to ~71 bp, p = 1.3 × 10⁻¹²⁹
over 202,078 CpG pairs.

**Numbers that moved.** The previously published region was
`chr6:28,633,534–28,633,601`, 6 CpGs, 67 bp, Šidák 1.2 × 10⁻¹⁰, −log₁₀ p 2.8–4.8,
Δβ −3.5 to −5.6 pp. It is superseded by:

| anchor | old | new |
|---|---|---|
| region | chr6:28,633,534–28,633,601 | chr6:28,633,493–28,633,701 |
| CpGs | 6 | 11 |
| span | 67 bp | 208 bp |
| Šidák *p* | 1.2 × 10⁻¹⁰ | 2.6 × 10⁻¹⁶ |
| `min_p` | — | 2.3 × 10⁻¹³ |
| −log₁₀ *p* range | 2.8 – 4.8 | 2.51 – 4.77 |
| Δβ range | −3.5 to −5.6 pp | −5.64 to −2.19 pp |
| best raw BACON *p* in region | 1.7 × 10⁻⁵ | 1.7 × 10⁻⁵ |

All eleven CpGs are hypomethylated in cases, all annotate to `Promoter`
(within ±1.5 kb of a TSS) and all to the same CpG-island north shore
(`N_Shore`) — cross-checked between the Zhou CGI knowledgebase
(`CGI.20220904.cm`, read with YAME) and the v41/v36 manifest, which agree. The
region's genes are `ENSG00000271440;ENSG00000287279`, distance to TSS −105 to
+101 bp. The summary bullet in chapter 08 was updated from the superseded
6-CpG / 1.2 × 10⁻¹⁰ claim.

The span is **208 bp**, computed as `max(CpG_beg) - min(CpG_beg) + 2` =
28,633,699 − 28,633,493 + 2. An earlier draft said 209; that came from reading
the region's exclusive end coordinate as inclusive. Corrected in the prose.

### 21.5 Narrative simplification

Per your decision, the DMR section is now a single default-seed run. The
"Seeing the machinery work" relaxed-seed subsection and the 1 × 10⁻⁸ demo run
were removed, along with the `dmr_demo/` output directory. The `min_p` callout
was kept — it corrects a real misreading, since `min_p` is the smallest
Stouffer–Liptak–Kechris-smoothed, FDR-corrected *regional* q among the region's
CpGs, not a raw probe p-value. The raw SLK p inside the region reaches
6.0 × 10⁻¹⁹, while the smallest raw BACON p among the eleven CpGs is
1.7 × 10⁻⁵ — three orders of magnitude apart, which is the point of the callout.

`08_dmr_locus.png` was regenerated against the 11-CpG region (9.0 × 5.2 in),
and `dmr_cpg_display.csv` / `dmr_cpg_annotation.csv` were regenerated from the
rebuilt `annotated.rds` rather than the stale CSV. The earlier version of the
table script had *inferred* annotations for the one CpG that v8.1 newly retains
and that the stale table therefore lacked; that inference is gone — all eleven
rows now carry real annotations from the manifest.

### 21.6 Stale shipped data files retired

Three files in `repo/data/` predated the v8.1 refit and were superseded rather
than regenerated, since nothing reads them:

- `06_ewas_bacon_toptable.csv` (uncompressed, Jul 30) — replaced by the
  `.csv.gz` the chapters actually `fread`
- `06_ewas_toptable.csv` (Jul 30) — pre-BACON, unreferenced
- `08_annotation/PTSD_ewas_annotated_results.csv.gz` (Jul 30)

All three were moved out of the repository tree. They are `.gitignore`d and were
never tracked, so no release tier or manifest entry changes.

### 21.7 Independent re-run of comb-p, and cleanup of stale run artefacts

The whole comb-p pipeline was re-run from the checked-in BED with exactly the
invocation the chapter prints, and every output compared byte-for-byte against
the shipped files:

```
comb-p pipeline -c 4 --seed 1e-4 --dist 200 --region-filter-p 0.05 \
    -p data/08_annotation/dmr/PTSD_dmr \
    data/08_annotation/PTSD_ewas_annotated_results.bed
```

`regions.bed.gz`, `regions-p.bed.gz`, `regions-t.bed`, `slk.bed.gz` and
`fdr.bed.gz` are all identical to what ships. The run also confirms two details
the chapter states implicitly: comb-p *derives* `--step 70` itself
(`calculated stepsize as: 70`, from the median inter-probe spacing), so the
printed command needs no `--step`; and the SLK step reports its own
`lambda: 1.02`, consistent with the chapter 06 BACON-adjusted λ = 1.023.

Two run-record files still described the superseded run and were replaced with
the current ones:

| file | was | now |
|---|---|---|
| `dmr/PTSD_dmr.args.txt` | `--seed 1e-15 … -p repo/data/… dmr_hdr.bed` | `--seed 1e-4 … -p data/08_annotation/dmr/PTSD_dmr data/08_annotation/PTSD_ewas_annotated_results.bed` |
| `dmr/PTSD_dmr.acf.txt` | correlation 0.05384, N = 202,081, p = 1.397 × 10⁻¹²⁹ | correlation 0.05384, N = 202,078, p = 1.349 × 10⁻¹²⁹ |

The ACF numbers moved because the old file came from `dmr_hdr.bed`, a scratch
BED with three more CpG pairs in range than the v8.1 BED the chapter now builds.
The correlation itself is unchanged to four figures; the chapter's rounded
"≈ 0.054 out to ~71 bp" text was already correct.

`dmr_cpg_annotation.csv` had been regenerated with a reduced column set
(`probe, chrm, beg, delta_beta, P.Value, bacon.p, feature, CGIposition`). It is
restored to its published schema — `probe, CpG_beg, CpG_end, genesUniq,
transcriptTypes, distTSS, feature, island, BIOS_eQTM_genes, delta_beta,
bacon.es, bacon.p` — now covering all eleven CpGs instead of seven. No number in
the chapter reads from this file; it is a supplementary download.

`dmr_source_comparison.csv` had picked up CRLF line endings; converted back to
LF. Its content had already been updated for the new coordinates
(chr6:28,633,493–28,633,701) and for the fact that COX8CP1 overlaps only the
region's first five CpGs rather than the whole region.

Removed from the repository: the superseded `dmr_demo/*` outputs, the duplicate
BED `PTSD_ewas_bed_from_v41.bed` (identical in construction to the shipped
`PTSD_ewas_annotated_results.bed`, which the chapter names), and two scratch
intermediates no chunk reads (`08_dist_data.rds`, `dmr_region_cpgs.rds`).
`regen_dmr_tables.R` no longer depends on a workspace scratch file
(`cgi_states.txt`); it takes the CGI relation from `annotated.rds` directly,
which is the same value the manifest carries.

## 22. Chapter 08 comb-p output documentation, `min_p` in the results table, and revised DMR prose

Author-revised prose for the comb-p section, plus the two additions that
revision called for. **No analysis was re-run and no published number moved.**
The values newly displayed were read out of the comb-p output already on disk.

### Added

- **A table of all seven files `comb-p pipeline` writes**, with the columns in
  each: `.args.txt` (invocation, version 0.50.6, run date), `.acf.txt`
  (`lag_min`, `lag_max`, `correlation`, `N`, `p`), `.slk.bed.gz` (`#chrom`,
  `start`, `end`, `p`, `region-p`), `.fdr.bed.gz` (adds `region-q`),
  `.regions.bed.gz` (headerless: chrom, start, end, `min_p`, `n_probes`),
  `.regions-p.bed.gz` (adds `z_p`, `z_sidak_p`), `.regions-t.bed` (the
  `--region-filter-p` survivors). Followed by a paragraph on which file to read.
- **`min_p` as a column in the candidate-regions table**, so the `min_p`
  callout that follows lands on a number the reader can see. New values shown,
  all read from `PTSD_dmr.regions-p.bed.gz`:
  - chr6:28,633,493–28,633,701 — `min_p` 2.3 × 10⁻¹³
  - chr3:138,608,544–138,608,573 — `min_p` 7.8 × 10⁻⁵
  - chr20:32,018,064–32,018,066 — `min_p` 2.1 × 10⁻⁵
- **A sentence stating the table is a subset**: five of seven columns, the three
  coordinate columns collapsed into one, and `z_p` omitted — with `z_p` for the
  chr6 region given as 3.6 × 10⁻²⁰.

### Changed (prose only)

- `min_p` callout: dropped the "This trips people up" opener.
- Seed callout: removed the 1 × 10⁻¹⁵-returns-zero-regions experiment; now
  directs the reader to justify whatever seed they choose.
- MHC callout: no longer credits the chapter 03 SNP mask with clearing the
  region; now names rare variants, structural variants, and undocumented
  genetic variation as live alternatives and points to causal-inference
  follow-up with genetic or expression data.
- ACF paragraph: dropped the "which is the premise the whole method rests on"
  clause.
- Closing paragraph: condensed from six lines to one sentence. The count reads
  **eleven** CpGs (the author's draft said "eight"; the region has eleven, as
  the figure, the tables, and `n_probes` in `regions-p.bed.gz` all agree).

### Verified, not changed — the 17 CpGs absent from the BED

Raised as a question about whether `end = CpG_beg + 2L` was costing probes. It
is not:

- All **756,251** tested probes have a Zhou v8.1 coordinate; none is missing.
- The **17** that do not reach the BED are removed by the primary-contig filter
  `grepl("^chr([0-9]+|X|Y)$", CpG_chrm)`: chr22_KI270879v1_alt 7,
  chr14_GL000009v2_random 3, chr1_KI270706v1_random 2, chr4_GL000008v2_random 2,
  chr14_KI270726v1_random 1, chr17_KI270857v1_alt 1, chr8_KI270821v1_alt 1.
- Best BACON p among the 17 is **0.037**; none reaches p < 1 × 10⁻⁵.
- `+2L` matches Zhou's own manifest convention — in `annotated.rds`,
  `CpG_end − CpG_beg` is 2 for 753,631 probes and 1 for 2,620.
- **A/B re-run of the full pipeline** confirms neither choice matters: with
  `end = beg + 1` (primary contigs) the chr6 region is 11 CpGs at
  `z_sidak_p` = 1.47 × 10⁻¹⁶; with `end = beg + 2` and *all* contigs included it
  is 11 CpGs at 2.59 × 10⁻¹⁶, byte-identical to what ships. Same three
  candidates in both. comb-p reads the end column only to reject spans
  > 100 kb (`cpv/_common.py:84`).

### Also confirmed — the figure-3 discrepancy

A local working copy showed a different Figure 3 than the site. The site is
correct: the local checkout was on commit `32ecf6d` (Aug 27) carrying the
**68,701-byte** superseded 6-CpG figure, while the deployed
`08_dmr_locus.png` is the **134,155-byte** 11-CpG rebuild. Resolved by pulling.
