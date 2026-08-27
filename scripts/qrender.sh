#!/usr/bin/env bash
# Quarto in the `methyl` conda env resolves its helper binaries under
# bin/tools/x86_64/, which is read-only here and holds only typst-gather.
# Every binary actually exists in bin/, so point quarto at them directly.
set -euo pipefail
E=/home/krferrier/.claude-science/conda/envs/methyl
W=/home/krferrier/.claude-science/orgs/5fc76295-8f5f-442d-a9d9-c1b932656a1a/workspaces/751c8b53-2493-44b0-acbe-1753dd0c4176

export PATH="$E/bin:$PATH"
export QUARTO_DENO="$E/bin/deno"
export QUARTO_PANDOC="$E/bin/pandoc"
export QUARTO_DART_SASS="$E/bin/sass"
export QUARTO_ESBUILD="$E/bin/esbuild"
export QUARTO_TYPST="$E/bin/typst"
export QUARTO_SHARE_PATH="$E/share/quarto"
export DENO_DOM_PLUGIN="$E/lib/deno_dom.so"
export DENO_INSTALL_ROOT="$E"
export R_LIBS="$W/repo/.r-libs/methyl:/home/krferrier/.claude-science/r-libs/751c8b53-2493-44b0-acbe-1753dd0c4176/methyl"

cd "$W/repo"
exec quarto "$@"
