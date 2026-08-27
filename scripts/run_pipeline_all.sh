set -e
W=/home/krferrier/.claude-science/orgs/5fc76295-8f5f-442d-a9d9-c1b932656a1a/workspaces/751c8b53-2493-44b0-acbe-1753dd0c4176
export R_LIBS="/home/krferrier/.claude-science/r-libs/751c8b53-2493-44b0-acbe-1753dd0c4176/methyl"
RS=/home/krferrier/.claude-science/conda/envs/methyl/bin/Rscript
cd $W/ewas_repo
D=$W/ewas_pipeline

echo "### pheno covariates in play"
head -1 $D/data/pheno.csv

echo "### combined EWAS (6 SVs + smoking proxy + array position, matching ch06)"
time $RS scripts/ewas.R --pheno $D/data/pheno.csv --methyl $D/data/mvals.csv.gz \
  --assoc PTSD --stratified no --chunk-size 5000 --processing-type multicore \
  --workers 8 --out-dir $D/run_grady/all --out-prefix all --out-type .csv.gz

echo "### stratify by sex"
$RS scripts/stratify.R --pheno $D/data/pheno.csv --methyl $D/data/mvals.csv.gz \
  --stratify sex --out-dir $D/run_grady/strat --threads 4

for S in F M; do
  echo "### ewas $S"
  time $RS scripts/ewas.R --pheno $D/run_grady/strat/$S/${S}_pheno.fst \
    --methyl $D/run_grady/strat/$S/${S}_mvals.fst --assoc PTSD --stratified yes \
    --chunk-size 5000 --processing-type multicore --workers 8 \
    --out-dir $D/run_grady/strat/$S --out-prefix $S --out-type .csv.gz
done

echo "### bacon x3"
$RS scripts/run_bacon.R -i $D/run_grady/all/PTSD_ewas_results.csv.gz \
  --out-dir $D/run_grady/all --out-type .csv.gz
for S in F M; do
  mkdir -p $D/run_grady/strat/$S/bacon_plots
  $RS scripts/run_bacon.R -i $D/run_grady/strat/$S/${S}_PTSD_ewas_results.csv.gz \
    --out-dir $D/run_grady/strat/$S --out-prefix $S --out-type .csv.gz
done

echo "### metal meta-analysis"
mkdir -p $D/run_grady/meta_analysis
bash scripts/metal_cmd.sh $D/run_grady/meta_analysis/PTSD_metal_commands.txt \
  $D/run_grady/meta_analysis/PTSD_ewas_meta_analysis_results_ \
  $D/run_grady/strat/F/F_PTSD_ewas_bacon_results.csv.gz \
  $D/run_grady/strat/M/M_PTSD_ewas_bacon_results.csv.gz
cd $D/run_grady/meta_analysis
$D/software/metal/build/metal/metal PTSD_metal_commands.txt | tail -20

echo "### outputs"
ls -la $D/run_grady/all/ $D/run_grady/strat/F/ $D/run_grady/strat/M/ $D/run_grady/meta_analysis/
