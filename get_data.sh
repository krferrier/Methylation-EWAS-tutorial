#!/usr/bin/env bash
# Fetch tutorial checkpoint data from Zenodo.
#
#   ./get_data.sh A_idats           # one tier
#   ./get_data.sh F_ewas_results H_annotation
#   ./get_data.sh all               # everything (~5.6 GB)
#
# A_idats holds the 192 raw IDAT files (96 samples x 2 channels) for the
# teaching subset, so chapter 01 can be run from raw data. Every other tier is
# a resume point that skips the steps before it.
#
# The record id below is the published one; override it by exporting ZENODO_RECORD
# (e.g. to pin an older version of the data).
#   Concept DOI (always newest): 10.5281/zenodo.22135215
#   This version:                10.5281/zenodo.22287946
set -euo pipefail

ZENODO_RECORD="${ZENODO_RECORD:-22287946}"
BASE="https://zenodo.org/records/${ZENODO_RECORD}/files"
ALL_TIERS=(A_idats B_qc C_normalized D_filtered E_model_inputs F_ewas_results G_pipeline_run H_annotation)

usage() { printf 'tiers: %s\n' "${ALL_TIERS[*]}"; exit 1; }
[ $# -ge 1 ] || usage
if [ "$1" = "all" ]; then want=("${ALL_TIERS[@]}"); else want=("$@"); fi

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
