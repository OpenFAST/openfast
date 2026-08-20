#!/bin/bash
# Prove the MirrorRotor guard rails actually fire (and only fire when they should).
# Uses the DEBUG build: these runs abort during init, nothing is baselined.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${OPENFAST_REPO:-$(cd "$HERE/../.." && pwd)}"
GC="${GUARDCHECK_DIR:-${TMPDIR:-/tmp}/openfast_guardcheck}"
mkdir -p "$GC"
RT=$REPO/reg_tests/r-test/glue-codes/openfast
OF=$REPO/build-docker-double-debug/glue-codes/openfast/openfast
SRC=$RT/5MW_Land_noDLL_Steady_MirrorRotor

# $1 = variant name, $2 = MirrorRotor (True/False), $3 = Linearize (True/False),
# $4 = Wake_Mod (1/3), $5 = CompAA (True/False)
mkvariant () {
  local name=$1 mirror=$2 lin=$3 wake=$4 aa=$5
  rm -rf "$GC/$name"; mkdir -p "$GC/$name"
  cp "$SRC"/*.fst "$SRC"/AeroDyn.dat "$SRC"/ElastoDyn.dat "$SRC"/ElastoDyn_Tower.dat "$GC/$name/"
  local fst="$GC/$name/5MW_Land_noDLL_Steady_MirrorRotor.fst"
  sed -i -E "s#^([[:space:]]*)[A-Za-z]+([[:space:]]+MirrorRotor)#\1$mirror\2#" "$fst"
  sed -i -E "s#^[A-Za-z]+([[:space:]]+Linearize)#$lin\1#" "$fst"
  sed -i -E "s#^[0-9]+([[:space:]]+Wake_Mod)#$wake\1#" "$GC/$name/AeroDyn.dat"
  sed -i -E "s#^[A-Za-z]+([[:space:]]+CompAA)#$aa\1#" "$GC/$name/AeroDyn.dat"
  if [ "$wake" = "3" ]; then
    cp "$RT/HelicalWake_OLAF"/*OLAF*.dat "$GC/$name/OLAF.dat" 2>/dev/null || \
      cp "$(find $RT/HelicalWake_OLAF -iname '*olaf*' -type f | head -1)" "$GC/$name/OLAF.dat"
    sed -i -E "s#^\"[^\"]*\"([[:space:]]+OLAFInputFileName)#\"OLAF.dat\"\1#" "$GC/$name/AeroDyn.dat"
  fi
  # keep it short - we only care about init
  sed -i -E "s#^([[:space:]]*)[0-9.]+([[:space:]]+TMax)#\1 0.5\2#" "$fst"
}

# $3 = "yes" if the message is expected, "no" if it must be absent
run () {
  local name=$1 expect=$2 want=$3
  cd "$GC/$name" || return
  timeout 300 "$OF" 5MW_Land_noDLL_Steady_MirrorRotor.fst > run.log 2>&1
  local rc=$?
  local got=no
  grep -qi "$expect" run.log && got=yes
  if [ "$got" = "$want" ]; then
    echo "PASS  $name  (rc=$rc)  guard fired=$got, expected=$want"
  else
    echo "FAIL  $name  (rc=$rc)  guard fired=$got, expected=$want"
    grep -i "FATAL\|Error" run.log | head -4 | sed 's/^/        /'
  fi
}

echo "=== guard rail proof, debug build ==="
mkvariant lin_mirror   True  True  1 False
mkvariant lin_cw       False True  1 False
mkvariant olaf_mirror  True  False 3 False
mkvariant olaf_cw      False False 3 False
mkvariant aa_mirror    True  False 1 True

run lin_mirror  "not yet supported with linearization"      yes
run lin_cw      "not yet supported with linearization"      no
run olaf_mirror "not yet supported with the OLAF"            yes
run olaf_cw     "not yet supported with the OLAF"            no
run aa_mirror   "not yet supported with the AeroAcoustics"   yes
