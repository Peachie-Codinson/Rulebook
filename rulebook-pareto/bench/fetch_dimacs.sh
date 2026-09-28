#!/usr/bin/env bash
# Fetch the DIMACS road networks used by the RA*pex paper.
#
# Only the first two objectives are public DIMACS files (distance and travel
# time). The third and fourth objectives, the query/instance files and the rules
# files live in the Google Drive folder linked from the reference repo's README
# and are NOT scriptable -- download those by hand into the same directory.
#
#   usage: bench/fetch_dimacs.sh [MAP] [OUTDIR]
#          MAP defaults to BAY (the map the paper reports on); NY is a smaller
#          one that is useful for a first smoke test.
set -euo pipefail

MAP="${1:-BAY}"
OUT="${2:-data/dimacs}"
BASE="${DIMACS_BASE:-https://www.diag.uniroma1.it/challenge9/data}"

mkdir -p "$OUT"
cd "$OUT"

# Verify the exact paths against the download page if these 404 -- the site has
# been reorganised before:
#   https://www.diag.uniroma1.it/challenge9/download.shtml
for kind in d t; do
    f="USA-road-${kind}.${MAP}.gr"
    if [ -f "$f" ]; then
        echo "have $f"
        continue
    fi
    echo "fetching $f.gz ..."
    curl -fSL --retry 3 -O "${BASE}/USA-road-${kind}/${f}.gz"
    gunzip -f "${f}.gz"
done

echo
echo "in $OUT:"
ls -la
cat <<EOF

Still needed, by hand, from the Google Drive folder linked in
https://github.com/Infus3d/Rulebook_approximation :

  USA-road-3.${MAP}.gr        third objective
  USA-road-4.${MAP}.gr        fourth objective (BAY only)
  ${MAP}_instances.txt        start/goal queries, one "s,t" per line
  <n>_rules_eps_<x>.txt       rulebook files

All .gr files must list the SAME arcs in the SAME order -- the loader checks
this arc for arc and tells you the first line where they diverge.
EOF
