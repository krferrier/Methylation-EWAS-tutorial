# Software environment and session information

The exact software state that produced every number, figure, and checkpoint in this
tutorial.

Generated: 2026-09-03

## Two ways to install: R alone, or conda

Everything except chapter 07 needs **only R**. `install_packages.R` in the repository
root installs the 4 CRAN and 18 Bioconductor packages the chapters use:

```bash
Rscript install_packages.R
```

It relies on `BiocManager`, which maps R 4.2 to Bioconductor 3.16 automatically. All 22
package names resolve in that release, and the script installs the pinned `matrixStats`
*before* anything else so a dependency resolution cannot bump it — `wateRmelon` and `sva`
both depend on `matrixStats` but neither declares a version floor, so `upgrade = "never"`
holds the pin.

Conda is required only for **chapter 07** (Snakemake) and the **comb-p section of chapter
08**. If you skip those, ignore the environment files entirely — chapter 06 already runs
the same association test directly in R.

## The three conda environments

If you are using conda, this tutorial uses three unique environments:

```bash
conda env create -f envs/methyl.yml     # chapters 01-06, 08 -- R 4.2.3
conda env create -f envs/smk.yml        # chapter 07         -- Snakemake 9.26.1
conda env create -f envs/combp.yml      # chapter 08 DMRs    -- Python 2.7.15
```

| Environment | Chapters | Core tool |
|---|---|---|
| `ewas-methyl` | 01-06, 08 | R 4.2.3, Bioconductor 3.16 |
| `ewas-smk` | 07 | Snakemake 9.26.1 |
| `ewas-combp` | 08 (DMR section) | comb-p on Python 2.7 |

Each file lists **only the packages the tutorial uses directly** and lets conda resolve
the dependency tree. Two version constraints are deliberate:

- **`r-base=4.2`** pins the Bioconductor generation to 3.16, which is what produced the
  checkpoints in the Zenodo record.
- **`r-matrixstats=1.0`** because `matrixStats` made `useNames = NA` defunct in 1.2.0
  (December 2023), while Bioconductor 3.16's `MatrixGenerics` still passes it. A newer
  `matrixStats` breaks `detectionP()` and `preprocessFunnorm()` with
  `useNames = NA is defunct`. `minfi` itself is not at fault — none of its own functions
  pass the argument.

Everything else floats, so the solver can pick builds that work on your platform.

## What is *not* in these files

**The Snakemake pipeline's R packages.** Chapter 07 runs `krferrier/EWAS`, which declares its own
per-rule environment (`envs/ewas.yaml`). Snakemake builds that itself on first run, because the workflow profile sets `software-deployment-method: conda`. 

**METAL.** The same pipeline compiles it from source in its `install_metal` rule.

**`EpiDISH`, `DMRcate`, `knowYourCG`.** These appear in the chapters as documented
alternatives inside non-executing (`eval: false`) blocks. Install them yourself if you
want to run those blocks.

## Troubleshooting: `pthread_create()` Error in chapter 02

If `preprocessFunnorm()` fails with:

```
ERROR; return code from pthread_create() is 22
```

`preprocessCore` maintains its own pthread pool, and the packaged build's threading path
fails on some multi-core Linux hosts and inside containers. It is **not universal** —
many people run the packaged build without trouble — so try the environment as shipped
first. If you hit it, rebuild without threading:

```bash
R CMD INSTALL --configure-args="--disable-threading" preprocessCore_1.60.2.tar.gz
```

The numbers here were produced with that rebuilt copy, which changes only the threading
strategy, not the arithmetic. Because it was installed into a user library rather than
via conda, it does not appear in `envs/methyl.yml`.

## Full `sessionInfo()`

Captured in `ewas-methyl` after attaching the 19 packages the chapters and
tutorial scripts execute. `sessionInfo()` reports **43 attached packages** rather
than 19, because attaching those pulls in their dependencies; a further 152
are loaded via namespace without being attached:

```
R version 4.2.3 (2023-03-15)
Platform: x86_64-conda-linux-gnu (64-bit)

Matrix products: default
BLAS/LAPACK: /home/krferrier/.claude-science/conda/envs/methyl/lib/libopenblasp-r0.3.33.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] preprocessCore_1.60.2                              
 [2] knitr_1.47                                         
 [3] QCEWAS_1.2-3                                       
 [4] sesame_1.16.0                                      
 [5] sesameData_1.16.0                                  
 [6] missMethyl_1.32.0                                  
 [7] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
 [8] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1 
 [9] methylGSA_1.16.0                                   
[10] bacon_1.26.0                                       
[11] ellipse_0.5.0                                      
[12] limma_3.54.0                                       
[13] sva_3.46.0                                         
[14] BiocParallel_1.32.5                                
[15] mgcv_1.9-1                                         
[16] nlme_3.1-165                                       
[17] genefilter_1.80.0                                  
[18] FlowSorted.Blood.EPIC_2.2.0                        
[19] ExperimentHub_2.6.0                                
[20] AnnotationHub_3.6.0                                
[21] BiocFileCache_2.6.0                                
[22] dbplyr_2.5.0                                       
[23] data.table_1.15.4                                  
[24] ggplot2_3.5.1                                      
[25] DT_0.33                                            
[26] minfi_1.44.0                                       
[27] bumphunter_1.40.0                                  
[28] locfit_1.5-9.9                                     
[29] iterators_1.0.14                                   
[30] foreach_1.5.2                                      
[31] Biostrings_2.66.0                                  
[32] XVector_0.38.0                                     
[33] SummarizedExperiment_1.28.0                        
[34] MatrixGenerics_1.10.0                              
[35] matrixStats_1.0.0                                  
[36] GenomicRanges_1.50.0                               
[37] GenomeInfoDb_1.34.9                                
[38] IRanges_2.32.0                                     
[39] S4Vectors_0.36.0                                   
[40] illuminaio_0.40.0                                  
[41] GEOquery_2.66.0                                    
[42] Biobase_2.58.0                                     
[43] BiocGenerics_0.44.0                                

loaded via a namespace (and not attached):
  [1] utf8_1.2.4                    tidyselect_1.2.1             
  [3] RSQLite_2.3.4                 AnnotationDbi_1.60.0         
  [5] htmlwidgets_1.6.4             grid_4.2.3                   
  [7] scatterpie_0.2.3              munsell_0.5.1                
  [9] codetools_0.2-20              statmod_1.5.0                
 [11] withr_3.0.0                   colorspace_2.1-0             
 [13] GOSemSim_2.24.0               filelock_1.0.3               
 [15] DOSE_3.24.0                   GenomeInfoDbData_1.2.9       
 [17] polyclip_1.10-6               farver_2.1.2                 
 [19] bit64_4.0.5                   rhdf5_2.42.0                 
 [21] downloader_0.4                treeio_1.22.0                
 [23] vctrs_0.6.5                   generics_0.1.3               
 [25] xfun_0.45                     gson_0.1.0                   
 [27] R6_2.5.1                      graphlayouts_1.1.0           
 [29] gridGraphics_0.5-1            bitops_1.0-7                 
 [31] rhdf5filters_1.10.0           cachem_1.1.0                 
 [33] reshape_0.8.9                 fgsea_1.24.0                 
 [35] DelayedArray_0.24.0           promises_1.3.0               
 [37] BiocIO_1.8.0                  scales_1.3.0                 
 [39] ggraph_2.1.0                  enrichplot_1.18.0            
 [41] gtable_0.3.5                  wheatmap_0.2.0               
 [43] tidygraph_1.3.0               rlang_1.1.4                  
 [45] splines_4.2.3                 lazyeval_0.2.2               
 [47] rtracklayer_1.58.0            BiocManager_1.30.23          
 [49] yaml_2.3.8                    reshape2_1.4.4               
 [51] GenomicFeatures_1.50.2        httpuv_1.6.15                
 [53] qvalue_2.30.0                 clusterProfiler_4.6.0        
 [55] tools_4.2.3                   ggplotify_0.1.2              
 [57] nor1mix_1.3-3                 RColorBrewer_1.1-3           
 [59] siggenes_1.72.0               Rcpp_1.0.12                  
 [61] plyr_1.8.9                    sparseMatrixStats_1.10.0     
 [63] progress_1.2.3                zlibbioc_1.44.0              
 [65] purrr_1.0.2                   RCurl_1.98-1.14              
 [67] prettyunits_1.2.0             openssl_2.2.0                
 [69] viridis_0.6.5                 cowplot_1.1.3                
 [71] ggrepel_0.9.5                 fs_1.6.4                     
 [73] magrittr_2.0.3                reactome.db_1.82.0           
 [75] patchwork_1.2.0               hms_1.1.3                    
 [77] mime_0.12                     xtable_1.8-4                 
 [79] HDO.db_0.99.1                 XML_3.99-0.17                
 [81] RobustRankAggreg_1.2.1        mclust_6.1                   
 [83] gridExtra_2.3                 compiler_4.2.3               
 [85] biomaRt_2.54.0                tibble_3.2.1                 
 [87] shadowtext_0.1.3              crayon_1.5.3                 
 [89] htmltools_0.5.8.1             ggfun_0.1.5                  
 [91] later_1.3.2                   tzdb_0.4.0                   
 [93] aplot_0.2.3                   tidyr_1.3.1                  
 [95] DBI_1.2.3                     tweenr_2.0.3                 
 [97] MASS_7.3-60.0.1               rappdirs_0.3.3               
 [99] Matrix_1.6-5                  readr_2.1.5                  
[101] cli_3.6.3                     quadprog_1.5-8               
[103] igraph_2.0.3                  pkgconfig_2.0.3              
[105] GenomicAlignments_1.34.0      xml2_1.3.6                   
[107] ggtree_3.6.0                  annotate_1.76.0              
[109] rngtools_1.5.2                multtest_2.54.0              
[111] beanplot_1.3.1                yulab.utils_0.1.4            
[113] doRNG_1.8.6                   scrime_1.3.5                 
[115] stringr_1.5.1                 digest_0.6.36                
[117] base64_2.0.1                  fastmatch_1.1-4              
[119] tidytree_0.4.6                edgeR_3.40.0                 
[121] DelayedMatrixStats_1.20.0     restfulr_0.0.15              
[123] curl_5.1.0                    shiny_1.8.1.1                
[125] Rsamtools_2.14.0              rjson_0.2.21                 
[127] jsonlite_1.8.8                lifecycle_1.0.4              
[129] Rhdf5lib_1.20.0               viridisLite_0.4.2            
[131] askpass_1.2.0                 fansi_1.0.6                  
[133] pillar_1.9.0                  lattice_0.22-6               
[135] KEGGREST_1.38.0               fastmap_1.2.0                
[137] httr_1.4.7                    survival_3.7-0               
[139] GO.db_3.16.0                  interactiveDisplayBase_1.36.0
[141] glue_1.7.0                    png_0.1-8                    
[143] BiocVersion_3.16.0            bit_4.0.5                    
[145] ggforce_0.4.2                 stringi_1.8.4                
[147] HDF5Array_1.26.0              blob_1.2.4                   
[149] org.Hs.eg.db_3.16.0           memoise_2.0.1                
[151] dplyr_1.1.4                   ape_5.8
```
