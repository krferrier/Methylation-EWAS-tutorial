#!/usr/bin/env bash
# Resume the v8.1 pipeline run from the BACON stage.
# The three EWAS arms (all / F / M) already completed; only bacon x3 and METAL remain.
set -euo pipefail
W=/home/krferrier/.claude-science/orgs/5fc76295-8f5f-442d-a9d9-c1b932656a1a/workspaces/751c8b53-2493-44b0-acbe-1753dd0c4176
D=$W/ewas_pipeline
RS=/home/krferrier/.claude-science/conda/envs/methyl/bin/Rscript
export R_LIBS=/home/krferrier/.claude-science/r-libs/751c8b53-2493-44b0-acbe-1753dd0c4176/methyl
unset R_LIBS_USER || true
cd "$W/ewas_repo"

echo "### bacon: combined arm"
time $RS scripts/run_bacon.R -i $D/run_grady/all/PTSD_ewas_results.csv.gz \
  --out-dir $D/run_grady/all --out-type .csv.gz

for S in F M; do
  echo "### bacon: $S arm"
  mkdir -p $D/run_grady/strat/$S/bacon_plots
  time $RS scripts/run_bacon.R -i $D/run_grady/strat/$S/${S}_PTSD_ewas_results.csv.gz \
    --out-dir $D/run_grady/strat/$S --out-prefix $S --out-type .csv.gz
done

echo "### metal meta-analysis"
mkdir -p $D/run_grady/meta_analysis
bash scripts/metal_cmd.sh $D/run_grady/meta_analysis/PTSD_metal_commands.txt \
  $D/run_grady/meta_analysis/PTSD_ewas_meta_analysis_results_ \
  $D/run_grady/strat/F/F_PTSD_ewas_bacon_results.csv.gz \
  $D/run_grady/strat/M/M_PTSD_ewas_bacon_results.csv.gz
cd $D/run_grady/meta_analysis
$D/software/metal/build/metal/metal PTSD_metal_commands.txt | tail -25

echo "### outputs"
ls -la $D/run_grady/all/ $D/run_grady/strat/F/ $D/run_grady/strat/M/ $D/run_grady/meta_analysis/
echo "### done"
