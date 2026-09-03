# Data manifest

Every file in the eight Zenodo tarballs, with its size and the tier's SHA-256 checksum. Paths inside the tarballs are relative to the repository root, with `repo/` rewritten to `tutorial/`, so extracting from the repo root puts each file where the `.qmd` documents expect it.

The published record is **https://doi.org/10.5281/zenodo.22135215** (concept DOI, always the newest
version; this build is record `22287946`, version DOI `10.5281/zenodo.22287946`). The
SHA-256 checksums below are what `get_data.sh` verifies after download; the record
page additionally lists Zenodo's own MD5 for each file.


### Tier `A_idats`

`ewas-tutorial-data-A_idats.tar.gz` — 1.30 GB compressed  
**Resume point:** ch01 — raw IDATs, run the tutorial from the beginning  
`sha256 322f789e56f3f6e0e5babd6316de77177f937dea9ae77c97a8a78acba3446769`

The 192 raw IDAT files for the 96-sample teaching subset (1,418,192,103 bytes uncompressed), two channels per sample, as downloaded from GEO. This is the only tier that is not a shortcut: every other tier lets you skip work, while this one is the raw input those shortcuts were computed from. Fetching it is equivalent to the per-sample download loops in [Setup](00_setup.qmd), and it is what lets `read.metharray()` in that chapter run.

First six files, of 192:

| file | bytes |
|---|---:|
| `repo/data/idats/GSM3853168_200932680028_R01C01_Grn.idat.gz` | 7,300,979 |
| `repo/data/idats/GSM3853168_200932680028_R01C01_Red.idat.gz` | 7,330,905 |
| `repo/data/idats/GSM3853169_200932680028_R02C01_Grn.idat.gz` | 7,320,887 |
| `repo/data/idats/GSM3853169_200932680028_R02C01_Red.idat.gz` | 7,336,603 |
| `repo/data/idats/GSM3853170_200932680028_R03C01_Grn.idat.gz` | 7,337,381 |
| `repo/data/idats/GSM3853170_200932680028_R03C01_Red.idat.gz` | 7,328,624 |
| … and 186 more | |

### Tier `B_qc`

`ewas-tutorial-data-B_qc.tar.gz` — 587 MB compressed  
**Resume point:** ch02 — skip reading 96 IDAT pairs (192 files)  
`sha256 e472da6aeb2dc188b90a2b4d6a184232cacb53c79c5f3bd846385da8a1528164`

| file | bytes |
|---|---:|
| `repo/data/01_RGset.rds` | 481,434,074 |
| `repo/data/01_detP.rds` | 108,980,124 |
| `repo/data/01_qc_pieces.rds` | 48,793 |
### Tier `C_normalized`

`ewas-tutorial-data-C_normalized.tar.gz` — 1,301 MB compressed  
**Resume point:** ch03 — skip funnorm normalization  
`sha256 e1f68b06b33cc7d74c6e5e4b0068b09998dce18f8f26b89189e6fc5bead777ae`

| file | bytes |
|---|---:|
| `repo/data/02_funnorm_grs.rds` | 1,252,923,324 |
| `repo/data/02_norm_pieces.rds` | 47,564,425 |
### Tier `D_filtered`

`ewas-tutorial-data-D_filtered.tar.gz` — 1053 MB compressed  
**Resume point:** ch04/ch05 — skip v8.1 mask filtering  
`sha256 220c7fbcaf2de0c758e5c9625140c101f3b50a10c3fd15dc7437106b6b9fd01e`

| file | bytes |
|---|---:|
| `repo/data/03_grs_filtered.rds` | 1,094,261,487 |
| `repo/data/03_mask_pieces.rds` | 3,691,143 |
| `repo/data/03_filter_funnel.rds` | 3,567,699 |
| `repo/data/EPIC.hg38.mask.v81.8code.tsv.gz` | 3,062,275 |

### Tier `E_model_inputs`

`ewas-tutorial-data-E_model_inputs.tar.gz` — 487 MB compressed  
**Resume point:** ch06 — skip cell composition, ComBat and SVA  
`sha256 0492d780b1ecd0e48a0b9b0b4bb75f3f664a03bc83274c848c9eb46c757b4b3c`

| file | bytes |
|---|---:|
| `repo/data/04_cc_full.rds` | 3,378 |
| `repo/data/04_validation.rds` | 8,114 |
| `repo/data/05_mvals_combat.rds` | 510,973,489 |
| `repo/data/05_sva.rds` | 11,687 |
| `repo/data/05_batch_pca.rds` | 14,241 |

### Tier `F_ewas_results`

`ewas-tutorial-data-F_ewas_results.tar.gz` — 112 MB compressed  
**Resume point:** ch07/ch08 — skip the limma + BACON run  
`sha256 54f2cc5e1241aef6f0ef2df175788c061922cba719c180e9ca8abcfa2b2ce9c3`

| file | bytes |
|---|---:|
| `repo/data/06_ewas.rds` | 56,035,031 |
| `repo/data/06_bacon_summary.rds` | 333 |
| `repo/data/06_ewas_bacon_toptable.csv.gz` | 61,029,526 |

### Tier `G_pipeline_run`

`ewas-tutorial-data-G_pipeline_run.tar.gz` — 265 MB compressed  
**Resume point:** ch07 — Snakemake pipeline outputs (combined + stratified + meta)  
`sha256 12eee2feff900d6630e271a8c62b8a8d4340d38a8f7b3ed0eb4d686f30bc62e5`

| file | bytes |
|---|---:|
| `ewas_pipeline/data/pheno.csv` | 17,485 |
| `ewas_pipeline/run_grady_all/PTSD_ewas_results.csv.gz` | 30,244,597 |
| `ewas_pipeline/run_grady_all/PTSD_ewas_bacon_results.csv.gz` | 57,673,185 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_fit.jpg` | 207,497 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_qqs.jpg` | 289,704 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_traces.jpg` | 584,236 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_posteriors.jpg` | 438,735 |
| `ewas_pipeline/run_grady/F/F_PTSD_ewas_results.csv.gz` | 29,286,839 |
| `ewas_pipeline/run_grady/F/F_PTSD_ewas_bacon_results.csv.gz` | 56,713,678 |
| `ewas_pipeline/run_grady/F/bacon_plots/F_PTSD_fit.jpg` | 193,550 |
| `ewas_pipeline/run_grady/F/bacon_plots/F_PTSD_qqs.jpg` | 266,265 |
| `ewas_pipeline/run_grady/F/bacon_plots/F_PTSD_traces.jpg` | 546,483 |
| … and 10 more | |

### Tier `H_annotation`

`ewas-tutorial-data-H_annotation.tar.gz` — 280 MB compressed  
**Resume point:** ch08 — Zhou manifests plus annotated results and comb-p outputs  
`sha256 3250d7c976ac69a6dcaa93d634ac575da72a7be5e58c0a8204462ae53f84191d`

| file | bytes |
|---|---:|
| `masks/EPIC.hg38.manifest.gencode.v41.tsv.gz` | 62,477,035 |
| `masks/EPIC.hg38.manifest.gencode.v36.tsv.gz` | 29,916,134 |
| `masks/EPIC.hg38.coord.tsv.gz` | 5,505,138 |
| `masks/EPIC.ordering.tsv.gz` | 7,756,709 |
| `repo/data/08_annotation/annotated.rds` | 69,949,757 |
| `repo/data/08_annotation/PTSD_ewas_annotated_results.bed` | 39,981,168 |
| `repo/data/08_annotation/PTSD_ewas_annotated_zhou.csv.gz` | 83,051,542 |
| `repo/data/08_annotation/bios_eqtm_hgnc_annotated.tsv` | 9,130,455 |
| `repo/data/08_annotation/08_dmr_combp_outputs.tar.gz` | 17,958,047 |

## Deliberately not distributed

These intermediates are written by the tutorial but never read back by any chapter. They are diagnostics whose summary CSVs and figures are committed to git.

| file | why it is excluded |
|---|---|
| `repo/data/05_mvals_combat_pos.rds` | redundant variant of 05_mvals_combat.rds (position-in-model diagnostic) |
| `repo/data/05_sv_worth_it.rds` | diagnostic sweep object; the CSV/PNG summaries ship in git |
| `repo/data/05_k_lambda_sweep.rds` | diagnostic sweep object; 05_k_lambda_sweep.csv ships in git |
| `repo/data/05_sva_prev.rds` | superseded SVA fit from an earlier design |
| `repo/data/05_position_comparison.rds` | diagnostic; summary CSV ships in git |
| `repo/data/05_k_sweep_posmod.rds` | diagnostic; summary CSV ships in git |
| `ewas_pipeline/data/mvals.csv.gz` | regenerable from tier D by scripts/build_pipeline_inputs.R (514 MB) |
