# Data manifest

Every file in the eight Zenodo tarballs, with its size and the tier's SHA-256 checksum. Paths inside the tarballs are relative to the repository root, with `repo/` rewritten to `tutorial/`, so extracting from the repo root puts each file where the `.qmd` documents expect it.

The published record is **https://doi.org/10.5281/zenodo.22135215** (concept DOI, always the newest
version; this build is record `23198791`, version DOI `10.5281/zenodo.23198791`). The
SHA-256 checksums below are what `get_data.sh` verifies after download; the record
page additionally lists Zenodo's own MD5 for each file.


### Tier `A_idats`

`ewas-tutorial-data-A_idats.tar.gz` — 1.30 GB compressed  
**Resume point:** chapter 1 — raw IDATs, run the tutorial from the beginning  
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
**Resume point:** chapter 2 — skip reading 96 IDAT pairs (192 files)  
`sha256 e472da6aeb2dc188b90a2b4d6a184232cacb53c79c5f3bd846385da8a1528164`

| file | bytes |
|---|---:|
| `repo/data/01_RGset.rds` | 481,434,074 |
| `repo/data/01_detP.rds` | 108,980,124 |
| `repo/data/01_qc_pieces.rds` | 48,793 |
### Tier `C_normalized`

`ewas-tutorial-data-C_normalized.tar.gz` — 1,301 MB compressed  
**Resume point:** chapter 3 (and chapter 5's EPISTRUCTURE example) — skip funnorm normalization  
`sha256 e1f68b06b33cc7d74c6e5e4b0068b09998dce18f8f26b89189e6fc5bead777ae`

| file | bytes |
|---|---:|
| `repo/data/02_funnorm_grs.rds` | 1,252,923,324 |
| `repo/data/02_norm_pieces.rds` | 47,564,425 |
### Tier `D_filtered`

`ewas-tutorial-data-D_filtered.tar.gz` — 1053 MB compressed  
**Resume point:** chapters 4 and 6 — skip v8.1 mask filtering  
`sha256 220c7fbcaf2de0c758e5c9625140c101f3b50a10c3fd15dc7437106b6b9fd01e`

| file | bytes |
|---|---:|
| `repo/data/03_grs_filtered.rds` | 1,094,261,487 |
| `repo/data/03_mask_pieces.rds` | 3,691,143 |
| `repo/data/03_filter_funnel.rds` | 3,567,699 |
| `repo/data/EPIC.hg38.mask.v81.8code.tsv.gz` | 3,062,275 |

`EPIC.hg38.mask.v81.8code.tsv.gz` is also tracked in git, so chapter 3 runs on a fresh
clone without fetching this tier. It stays listed here because the deposited tarball
contains it; if you fetch the tier you simply overwrite the file with an identical copy.

### Tier `E_model_inputs`

`ewas-tutorial-data-E_model_inputs.tar.gz` — 482 MB compressed  
**Resume point:** chapter 7 — skip cell composition, ComBat and SVA  
`sha256 a4ce8854c8167db5a7f1f99dd1201b285127c26715907d2c4bdc0fc71c45e771`

| file | bytes |
|---|---:|
| `repo/data/04_cc_full.rds` | 3,378 |
| `repo/data/04_validation.rds` | 8,203 |
| `repo/data/06_mvals_combat.rds` | 505,149,587 |
| `repo/data/06_sva.rds` | 14,422 |
| `repo/data/06_batch_pca.rds` | 18,339 |

### Tier `F_ewas_results`

`ewas-tutorial-data-F_ewas_results.tar.gz` — 100 MB compressed  
**Resume point:** chapters 8 and 9 — skip the limma + BACON run  
`sha256 b1660b8ea52048adfb27d6f3acaaec943ca09b2cc687140e800e8999ccd41468`

| file | bytes |
|---|---:|
| `repo/data/07_ewas.rds` | 50,473,587 |
| `repo/data/07_bacon_summary.rds` | 187 |
| `repo/data/07_ewas_bacon_toptable.csv.gz` | 53,942,729 |

### Tier `G_pipeline_run`

`ewas-tutorial-data-G_pipeline_run.tar.gz` — 265 MB compressed  
**Resume point:** chapter 8 — Snakemake pipeline outputs (combined + stratified + meta)  
`sha256 c67fcbe9e7a9cfaf342e710f790d79ba4ea8bbca95cda3510c9a8cee0ef7be85`

| file | bytes |
|---|---:|
| `ewas_pipeline/data/pheno.csv` | 20,219 |
| `ewas_pipeline/run_grady_all/PTSD_ewas_results.csv.gz` | 30,244,197 |
| `ewas_pipeline/run_grady_all/PTSD_ewas_bacon_results.csv.gz` | 57,795,231 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_fit.jpg` | 201,552 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_qqs.jpg` | 286,681 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_traces.jpg` | 563,129 |
| `ewas_pipeline/run_grady_all/bacon_plots/PTSD_posteriors.jpg` | 443,820 |
| `ewas_pipeline/run_grady/F/F_PTSD_ewas_results.csv.gz` | 29,265,907 |
| `ewas_pipeline/run_grady/F/F_PTSD_ewas_bacon_results.csv.gz` | 56,487,435 |
| `ewas_pipeline/run_grady/F/bacon_plots/F_PTSD_fit.jpg` | 196,995 |
| `ewas_pipeline/run_grady/F/bacon_plots/F_PTSD_qqs.jpg` | 264,950 |
| `ewas_pipeline/run_grady/F/bacon_plots/F_PTSD_traces.jpg` | 544,158 |
| … and 10 more | |

### Tier `H_annotation`

`ewas-tutorial-data-H_annotation.tar.gz` — 177 MB compressed  
**Resume point:** chapter 9 — annotated results and comb-p outputs  
`sha256 ece3e940f87d12e5946d40f673cba1a8e2f8b355f7bc1f154c3863b463c77385`

| file | bytes |
|---|---:|
| `repo/data/09_annotation/annotated.rds` | 76,578,899 |
| `repo/data/09_annotation/PTSD_ewas_annotated_results.bed` | 39,982,327 |
| `repo/data/09_annotation/PTSD_ewas_annotated_zhou.csv.gz` | 75,766,376 |
| `repo/data/09_annotation/dmr/PTSD_dmr.slk.bed.gz` | 8,362,568 |
| `repo/data/09_annotation/dmr/PTSD_dmr.fdr.bed.gz` | 9,573,592 |

## Deliberately not distributed

| file | why it is excluded |
|---|---|
| `ewas_pipeline/data/mvals.csv.gz` | regenerable from tier E by chapter 8's pipeline-input code (514 MB) |
