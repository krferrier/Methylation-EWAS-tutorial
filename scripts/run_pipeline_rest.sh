set -e
W=/home/krferrier/.claude-science/orgs/5fc76295-8f5f-442d-a9d9-c1b932656a1a/workspaces/751c8b53-2493-44b0-acbe-1753dd0c4176
export R_LIBS="$W/repo/.r-libs/methyl:/home/krferrier/.claude-science/r-libs/751c8b53-2493-44b0-acbe-1753dd0c4176/methyl"
RS=/home/krferrier/.claude-science/conda/envs/methyl/bin/Rscript
cd $W/ewas_repo
D=$W/ewas_pipeline

echo "### stratify"
$RS scripts/stratify.R --pheno $D/data/pheno.csv --methyl $D/data/mvals.csv.gz \
  --stratify sex --out-dir $D/run_grady/strat --threads 4
ls -la $D/run_grady/strat/*/

for S in F M; do
  echo "### ewas $S"
  $RS scripts/ewas.R --pheno $D/run_grady/strat/$S/${S}_pheno.fst \
    --methyl $D/run_grady/strat/$S/${S}_mvals.fst --assoc PTSD --stratified yes \
    --chunk-size 5000 --processing-type multicore --workers 8 \
    --out-dir $D/run_grady/strat/$S --out-prefix $S --out-type .csv.gz
done

for R in all strat/F strat/M; do
  echo "### bacon $R"
  $RS scripts/run_bacon.R -i $D/run_grady/$R/*_PTSD_ewas_results.csv.gz \
    --out-dir $D/run_grady/$R --out-type .csv.gz
done

echo "### metal"
bash scripts/metal_cmd.sh $D/run_grady/meta_analysis/PTSD_metal_commands.txt \
  $D/run_grady/meta_analysis/PTSD_ewas_meta_analysis_results_ \
  $D/run_grady/strat/F/F_PTSD_ewas_bacon_results.csv.gz \
  $D/run_grady/strat/M/M_PTSD_ewas_bacon_results.csv.gz
cd $D/run_grady/meta_analysis && $W/ewas_pipeline/software/metal/build/metal/metal PTSD_metal_commands.txt | tail -25
ls -la $D/run_grady/meta_analysis/
