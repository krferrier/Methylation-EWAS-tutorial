# Data manifest

Every file in the seven Zenodo tarballs, with its size and the tier's SHA-256 checksum. Paths inside the tarballs are relative to the repository root, with `repo/` rewritten to `tutorial/`, so extracting from the repo root puts each file where the `.qmd` documents expect it.

### Tier `B_qc`

`ewas-tutorial-data-B_qc.tar.gz` — 587 MB compressed  
**Resume point:** ch02 — skip reading 622 IDAT pairs  
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

`ewas-tutorial-data-D_filtered.tar.gz` — 1,102 MB compressed  
**Resume point:** ch04/ch05 — skip v8.1 mask filtering  
`sha256 174ac17ee79f50b9efa55231f982a63d529765e89e8e0e2ff94b492764e359c9`

| file | bytes |
|---|---:|
| `repo/data/03_grs_filtered.rds` | 1,094,261,487 |
| `repo/data/03_mask_pieces.rds` | 3,691,143 |
| `repo/data/03_filter_funnel.rds` | 3,567,699 |
### Tier `E_model_inputs`

`ewas-tutorial-data-E_model_inputs.tar.gz` — 511 MB compressed  
**Resume point:** ch06 — skip cell composition, ComBat and SVA  
`sha256 989a7bfb6d161f107cac7c9cd3aaab754acb88d49e923a8d426d1d352137506d`

| file | bytes |
|---|---:|
| `repo/data/04_cc_full.rds` | 3,378 |
| `repo/data/04_validation.rds` | 8,114 |
| `repo/data/05_mvals_combat.rds` | 510,973,400 |
| `repo/data/05_sva.rds` | 11,687 |
| `repo/data/05_batch_pca.rds` | 14,241 |
### Tier `F_ewas_results`

`ewas-tutorial-data-F_ewas_results.tar.gz` — 117 MB compressed  
**Resume point:** ch07/ch08 — skip the limma + BACON run  
`sha256 a578bff697071d6bf099cf0410934e3d9432776f0ac76ef82feb18bbacb04845`

| file | bytes |
|---|---:|
| `repo/data/06_ewas.rds` | 56,036,353 |
| `repo/data/06_bacon_summary.rds` | 333 |
| `repo/data/06_ewas_bacon_toptable.csv.gz` | 60,794,109 |
### Tier `G_pipeline_run`

`ewas-tutorial-data-G_pipeline_run.tar.gz` — 220 MB compressed  
**Resume point:** ch07 — Snakemake pipeline outputs (no need to re-run 4 h of EWAS)  
`sha256 5de981a73c0126e41be33bfe37a24b2c319bf96f4ad734e2a979b93cc6406374`

| file | bytes |
|---|---:|
| `ewas_pipeline/data/pheno.csv` | 17,485 |
| `ewas_pipeline/run_grady/all/PTSD_ewas_results.csv.gz` | 30,245,000 |
| `ewas_pipeline/run_grady/all/PTSD_ewas_bacon_results.csv.gz` | 57,674,325 |
| `ewas_pipeline/run_grady/strat/F/F_PTSD_ewas_bacon_results.csv.gz` | 56,714,625 |
| `ewas_pipeline/run_grady/strat/M/M_PTSD_ewas_bacon_results.csv.gz` | 56,700,202 |
| `ewas_pipeline/run_grady/meta_analysis/PTSD_ewas_meta_analysis_results_1.txt` | 42,987,000 |
| `ewas_pipeline/run_grady/meta_analysis/PTSD_metal_commands.txt` | 774 |
### Tier `H_annotation`

`ewas-tutorial-data-H_annotation.tar.gz` — 305 MB compressed  
**Resume point:** ch08 — Zhou manifests plus annotated results and comb-p output  
`sha256 d82cda05e7f4be344496c778614bab09f6a5b5da70dcffe13c4a1ef9b0ca218b`

| file | bytes |
|---|---:|
| `masks/EPIC.hg38.manifest.gencode.v41.tsv.gz` | 62,477,035 |
| `masks/EPIC.hg38.manifest.gencode.v36.tsv.gz` | 29,916,134 |
| `masks/EPIC.hg38.coord.tsv.gz` | 5,505,138 |
| `repo/data/08_annotation/annotated.rds` | 69,948,881 |
| `repo/data/08_annotation/PTSD_ewas_annotated_results.bed` | 39,982,171 |
| `repo/data/08_annotation/PTSD_ewas_annotated_zhou.csv.gz` | 83,052,368 |
| `repo/data/08_annotation/bios_eqtm_hgnc_annotated.tsv` | 9,130,455 |
| `repo/data/08_annotation/08_dmr_combp_outputs.tar.gz` | 34,346,704 |

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
