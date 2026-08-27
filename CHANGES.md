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
