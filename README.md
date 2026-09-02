# An EWAS walkthrough: from raw IDATs to functional annotation

A reproducible, eleven-chapter Quarto tutorial that takes a publicly available Illumina
Infinium MethylationEPIC dataset from raw IDAT files through quality control,
normalization, probe filtering, cell-composition estimation, batch-effect diagnosis,
an epigenome-wide association study, a Snakemake production pipeline, and functional
annotation of the results.

**Rendered tutorial:** https://krferrier.github.io/Methylation-EWAS-tutorial/

**Data:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22135215.svg)](https://doi.org/10.5281/zenodo.22135215)

## The dataset

Grady Trauma Project, NCBI GEO accession **GSE132203** — whole-blood EPIC methylation
with matched phenotype data. The walkthrough uses a curated subset of **96 arrays**
(approximately 95% African American). All 96 pass QC and are carried through
normalization, probe filtering, cell-composition estimation, and batch diagnosis;
the association model is fit on the **87** samples with recorded PTSD status
(32 cases, 55 controls) — nine samples have no PTSD status in the deposit and drop
out when the exposure enters the model, which is a problem of missingness of the phenotype of interest rather than a
quality failure. After filtering under the Zhou EPIC `mask.cm` **v8.1** recommended
mask, **756,251 CpGs** are tested.

The point of the walkthrough is the method and its calibration (λ ≈ 1.02), not a
discovery. Treat the association results as a teaching artifact, not a finding.

## Chapters

| | chapter | what it covers                                                                 |
|---|---|--------------------------------------------------------------------------------|
| 00 | `00_setup.qmd` | environment, packages, CpG-island vocabulary                                   |
| 00b | `00b_dataset_catalog.qmd` | how the GEO dataset was chosen from 70 candidates                              |
| 01 | `01_qc.qmd` | IDAT import, detection *p*, sample and probe QC                                |
| 02 | `02_normalization.qmd` | functional normalization                                                       |
| 03 | `03_probe_filtering.qmd` | masking, detection *p*, sex chromosomes                                          |
| 04 | `04_cell_composition.qmd` | IDOL reference-based deconvolution                                             |
| 05 | `05_batch_effects.qmd` | PCA, sex-stratified ComBat, SVA (k = 6)                                        |
| 06 | `06_ewas.qmd` | limma EWAS, BACON bias and inflation correction                                |
| 07 | `07_pipeline.qmd` | the same analysis as a Snakemake pipeline, sex-stratified, METAL meta-analysis |
| 08 | `08_annotation.qmd` | gene and CGI annotation, gometh / methylGSA / KYCG enrichment, comb-p DMRs     |

## Repository layout

```
tutorial/          the .qmd chapters, _quarto.yml, _setup.R, references.bib
tutorial/data/     small committed data: figures, summary CSVs, sample sheets
scripts/           the R and shell scripts used to (re)compute each checkpoint
docs/              the rendered site, served by GitHub Pages
get_data.sh        fetches the large checkpoints from Zenodo
MANIFEST.md        every distributed file, its size, and per-tier checksums
dev-notes/         maintainer notes: analysis decisions, change log, release checklist
                   (not needed to read, render, or reuse the tutorial)
```

## Where the data is stored

The tutorial tests the full post-QC CpG set, so the
intermediate objects are large: the normalized `GenomicRatioSet` alone is 1.2 GiB and the
whole checkpoint set is about 4.3 GiB on disk. So, **code, prose and every file under 5 MB are in git; the large checkpoints are a
versioned Zenodo record with a DOI.** Zenodo allows 50 GB per record, is free, has no
bandwidth metering, and gives the dataset a citable identifier — which is the right
outcome for a teaching resource anyway.

## Getting the data

You do not need all of it. Each tarball is a **resume point**: download only the tier for
the chapter you want to start from, and every later chapter will run.

| tier | download | resume at | contents |
|---|---:|---|---|
| `B_qc` | 587 MB | **ch02** | `RGChannelSet`, detection *p* matrix, QC summaries |
| `C_normalized` | 1301 MB | **ch03** | funnorm `GenomicRatioSet` |
| `D_filtered` | 1102 MB | **ch04** | mask-filtered `GenomicRatioSet` (756,273 probes), mask pieces, filter funnel |
| `E_model_inputs` | 511 MB | **ch06** | cell proportions, ComBat M-values, SVA fit (k = 6), batch PCA |
| `F_ewas_results` | 117 MB | **ch07** | limma top table, BACON summary, corrected top table |
| `G_pipeline_run` | 220 MB | **ch07 figures** | Snakemake pipeline outputs: combined arm, F/M strata, METAL meta-analysis |
| `H_annotation` | 306 MB | **ch08** | Zhou gencode v36/v41 manifests, annotated results, comb-p output |

```bash
git clone https://github.com/krferrier/Methylation-EWAS-tutorial.git
cd Methylation-EWAS-tutorial

./get_data.sh F_ewas_results H_annotation   # enough to run ch07 and ch08
./get_data.sh all                           # everything, ~4.1 GB
```

`get_data.sh` verifies each tarball against `SHA256SUMS.txt` before extracting, and
extracts in place so files land where the `.qmd` documents look for them.

Starting from chapter 01 instead means downloading the 622 raw IDAT files from GEO
yourself; `01_qc.qmd` documents how.

## Rendering

```bash
bash qrender.sh render                       # whole site -> tutorial/_site
bash qrender.sh render 06_ewas.qmd           # one chapter
```

Chapter 04 needs `FlowSorted.Blood.EPIC`, which is a large Bioconductor experiment
package; if it is not installed, use tier `E_model_inputs` and the chapter will read the
cached proportions instead of recomputing them.

## Reusing the pipeline

Chapter 07 drives the Snakemake EWAS pipeline at
**https://github.com/krferrier/EWAS**. Clone it as a sibling directory
(`../ewas_pipeline`) if you want to run it yourself; tier `G_pipeline_run` contains its
outputs so the chapter renders without recomputing.

## Citing

**If you use this tutorial — in teaching, in a rotation project, in a methods
section, or as the basis for your own analysis — please cite it.** The MIT
license does not require this, but academic credit is how work like this is
sustained, and it costs you one line.

There are two things to cite, and which you need depends on what you used:

- **The tutorial itself** (prose, code, or the approach) — cite the repository.
- **The checkpoint data** (anything fetched by `get_data.sh`) — cite the Zenodo
  record as well.

### Cite the tutorial

> Ferrier, K. (2026). *An EWAS walkthrough: from raw IDATs to functional
> annotation.* https://github.com/krferrier/Methylation-EWAS-tutorial

BibTeX:

```bibtex
@misc{ferrier2026ewaswalkthrough,
  author       = {Ferrier, Kendra},
  title        = {An {EWAS} walkthrough: from raw {IDAT}s to functional annotation},
  year         = {2026},
  howpublished = {\url{https://github.com/krferrier/Methylation-EWAS-tutorial}},
  note         = {Rendered at \url{https://krferrier.github.io/Methylation-EWAS-tutorial/}}
}
```

GitHub also reads `CITATION.cff` in this repository — use the **"Cite this
repository"** button in the sidebar to get APA or BibTeX directly.

### Cite the checkpoint data

> Ferrier, K. (2026). *EWAS tutorial (GSE132203, Zhou EPIC mask.cm v8.1):
> checkpoint data* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.22135215

```bibtex
@dataset{ferrier2026ewasdata,
  author    = {Ferrier, Kendra},
  title     = {{EWAS} tutorial ({GSE132203}, {Zhou EPIC} mask.cm v8.1): checkpoint data},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22135215}
}
```

That is the *concept* DOI: it always resolves to the newest version of the data.
To cite the exact files this tutorial was built from, use the version DOI
https://doi.org/10.5281/zenodo.22135216 (record `22135216`).

### Also cite the upstream sources

The tutorial is a walkthrough over other people's data and methods. If your work
depends on them, cite them directly rather than by way of this tutorial:

- **The dataset** — Grady Trauma Project, NCBI GEO
  [GSE132203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132203). The
  methylation data belongs to its original depositors; it appears here only as
  derived intermediates for teaching.
- **The probe mask** — Zhou, W., Laird, P. W., & Shen, H. (2017). Comprehensive
  characterization, annotation and innovative use of Infinium DNA methylation
  BeadChip probes. *Nucleic Acids Research*, 45(4), e22.
  https://doi.org/10.1093/nar/gkw967
- **The analysis tools** — `minfi`, `sesame`, `limma`, `sva`, `bacon`, `comb-p`,
  `missMethyl`, `METAL` and the rest are cited in each chapter's reference list.

## License

**MIT** (see [`LICENSE`](LICENSE)) — for the tutorial's code and prose. You may
use, modify, and redistribute it, including commercially, provided the copyright
notice is preserved. Citation is requested as a scholarly courtesy (see
[Citing](#citing)), not required by the license.

Two things the MIT license here does **not** cover:

- **The GSE132203 methylation data** is the property of its original depositors
  and is redistributed only as derived intermediates for teaching. Its use is
  governed by the terms of the GEO deposit, not by this repository's license.
- **The software the tutorial calls** — Bioconductor packages, comb-p, METAL and
  so on — carries each project's own license.
