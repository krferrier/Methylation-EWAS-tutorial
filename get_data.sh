#!/usr/bin/env bash
# Fetch tutorial checkpoint data from Zenodo.
#
#   ./get_data.sh B_qc              # one tier
#   ./get_data.sh F_ewas_results H_annotation
#   ./get_data.sh all               # everything (~4.1 GB)
#
# Set ZENODO_DOI once the record is published, or export it in your shell.
set -euo pipefail

ZENODO_RECORD="${ZENODO_RECORD:-REPLACE_WITH_RECORD_ID}"
BASE="https://zenodo.org/records/${ZENODO_RECORD}/files"
ALL_TIERS=(B_qc C_normalized D_filtered E_model_inputs F_ewas_results G_pipeline_run H_annotation)

usage() { printf 'tiers: %s\n' "${ALL_TIERS[*]}"; exit 1; }
[ $# -ge 1 ] || usage
if [ "$1" = "all" ]; then want=("${ALL_TIERS[@]}"); else want=("$@"); fi

if [ "$ZENODO_RECORD" = "REPLACE_WITH_RECORD_ID" ]; then
  echo "ERROR: set ZENODO_RECORD to the numeric Zenodo record id first." >&2; exit 2
fi

for t in "${want[@]}"; do
  case " ${ALL_TIERS[*]} " in *" $t "*) ;; *) echo "unknown tier: $t" >&2; usage ;; esac
  tar="ewas-tutorial-data-${t}.tar.gz"
  if [ ! -f "$tar" ]; then
    echo "--> downloading $tar"
    curl -fL --retry 3 -o "$tar" "${BASE}/${tar}?download=1"
  fi
  echo "--> verifying $tar"
  grep " ${tar}\$" SHA256SUMS.txt | sha256sum -c -
  echo "--> extracting $tar"
  tar -xzf "$tar"
  echo "    ok: $t"
done
echo "Done. Data landed under tutorial/data/ (and ewas_pipeline/, masks/ for tiers G and H)."
