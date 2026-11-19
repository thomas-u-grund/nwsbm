#!/bin/zsh
# nwsbm test runner (macOS). Usage: cscripts/run_tests.sh [jobs]. Run from anywhere; the repository root is the parent of cscripts/.
# NWCORE: directory of the core nwcommands package (default: the owner's checkout); NWSBM_SCRATCH: scratch directory.
# Each test is an independent do-file run in its own Stata process (-e, never -b); logs are test_<name>.log in the repository root.
cd "$(dirname "$0")/.." || exit 1
python3 tools/check_version.py || exit 1
jobs=${1:-3}
STATA=${STATA:-"/Applications/Stata 17.0/StataBE.app/Contents/MacOS/StataBE"}
tests=($(cd cscripts; ls test_*.do | sed "s/\.do$//"))
running=0
for t in $tests; do
  ( NWSBM_TAG=$t "$STATA" -e do cscripts/$t.do > /dev/null 2>&1
    printf '%-28s PASS lines: %s   errors: %s\n' $t "$(grep -c -E '^PASS' $t.log)" "$(grep -E '^r\([0-9]+\)' $t.log | sort -u | tr '\n' ' ')" ) &
  running=$((running + 1))
  if (( running >= jobs )); then wait -n 2>/dev/null || wait; running=$((running - 1)); fi
done
wait
echo "done (logs: test_*.log in the repository root; a test without a log counts as NOT RUN)"
