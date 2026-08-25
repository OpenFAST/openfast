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

# The deck refers to ../5MW_Baseline for the inflow and airfoil data.  Without this
# the runs abort while reading their inputs, before reaching the guard at all - and a
# guard that is never reached looks exactly like a guard that stayed silent, so the
# clockwise controls would pass for the wrong reason.
ln -sfn "$RT/5MW_Baseline" "$GC/5MW_Baseline"

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
run aa_mirror   "not yet supported with the AeroAcoustics"   yes

# The free wake is supported, so these two check that the restriction stays gone
# rather than that it fires. If someone reintroduces a refusal for OLAF, both fail.
run olaf_mirror "not yet supported with the OLAF"            no
run olaf_cw     "not yet supported with the OLAF"            no

# The VTK blade surfaces are built from the airfoil coordinate files rather than
# from the mesh, so they carry a sign of their own that no channel comparison and
# no node comparison can reach. This uses the release build and takes about 12 s.
echo
echo "=== VTK surface mirror ==="
python3 "$HERE/check_vtk_surface_mirror.py"

# ---------------------------------------------------------------------------
# FAST.Farm guards
#
# FAST.Farm is a separate executable and "make openfast" does not rebuild it.
# That is not hypothetical: when this check was written the binary in
# build-docker-double was six days stale and still carried a guard message that
# had been removed from the branch, which very nearly recorded a restriction
# that no longer exists.  A check that runs FAST.Farm catches that.
#
# The debug build does not include FAST.Farm, so fall back to the release one.
FFBIN=""
for cand in "$REPO/build-docker-double-debug/glue-codes/fast-farm/FAST.Farm" \
            "$REPO/build-docker-double/glue-codes/fast-farm/FAST.Farm"; do
  [ -x "$cand" ] && FFBIN="$cand" && break
done

# Prefer the staged build tree: the DISCON controller libraries are built, not
# stored, so they do not exist in the r-test source tree at all, and a deck that
# cannot load its controller aborts before reaching any guard.
FRT=$REPO/reg_tests/r-test/glue-codes/fast-farm
FSTAGE=$REPO/build-docker-double/reg_tests/glue-codes/fast-farm
[ -d "$FSTAGE/5MW_Baseline/ServoData" ] && FSRC=$FSTAGE || FSRC=$FRT
FG=$GC/farm

# $1 = variant name.  Edits are applied to turbine 1 only; turbine 2 is the
# control within every run, so a guard that fires for the wrong turbine shows up.
mkfarm () {
  local name=$1
  rm -rf "$FG/$name"; mkdir -p "$FG"
  cp -r "$FSRC/TSinflow" "$FG/$name"
  # Long enough to get through initialisation, short enough not to be a test run.
  sed -i -E 's/^[0-9.]+([[:space:]]+TMax)/5.0\1/' "$FG/$name/FAST.Farm.fstf"
}

# $1 = variant, $2 = text to look for, $3 = "yes" if it must appear
runfarm () {
  local name=$1 expect=$2 want=$3
  if [ -z "$FFBIN" ]; then
    echo "SKIP  $name  (no FAST.Farm binary; build it with 'make FAST.Farm')"
    return
  fi
  cd "$FG/$name" || return
  timeout 900 "$FFBIN" FAST.Farm.fstf > run.log 2>&1
  local rc=$? got=no
  grep -qF "$expect" run.log && got=yes
  if [ "$got" = "$want" ]; then
    echo "PASS  $name  (rc=$rc)  guard fired=$got, expected=$want"
  else
    echo "FAIL  $name  (rc=$rc)  guard fired=$got, expected=$want"
    grep -iE "FATAL|Error" run.log | head -4 | sed 's/^/        /'
  fi
}

echo
echo "=== FAST.Farm guard rails ==="

# The siblings the decks reach for.  Without them the runs abort while reading
# inputs, before reaching any guard, and a control that never initialises looks
# exactly like a control whose guard stayed silent.
mkdir -p "$FG"
ln -sfn "$FSRC/WAT_MannBoxDB" "$FG/WAT_MannBoxDB"
ln -sfn "$FSRC/5MW_Baseline"  "$FG/5MW_Baseline"

MIRROR_MSG="MirrorRotor is not yet supported with FAST.Farm."
NROTOR_MSG="Only one rotor per OpenFAST instance is supported with FAST.Farm."

# Control: one rotor, not mirrored.  Neither guard may fire, and it must run.
mkfarm farm_control
runfarm farm_control "$MIRROR_MSG" no
runfarm farm_control "$NROTOR_MSG" no
runfarm farm_control "FAST.Farm terminated normally." yes

# Multiple rotors in one instance.  Every rotor reference in FASTWrapper is
# rotors(1), so rotors 2 and beyond would be simulated as though absent.
mkfarm farm_multirotor
sed -i -E 's/^([[:space:]]*)1([[:space:]]+NRotors)/\1          2\2/' "$FG/farm_multirotor/FFTest_WT1.fst"
sed -i -E 's/^([[:space:]]*)F([[:space:]]+MirrorRotor)/\1      F F\2/' "$FG/farm_multirotor/FFTest_WT1.fst"
# NRotors > 1 needs a second block of per-rotor input files or the deck will not
# read, and a deck that fails to read never reaches the guard.
awk '/^---------------------- OUTPUT/ && !ins {
       print "---------------------- INPUT FILES Rotor 2 -------------------------------------";
       print "\"NRELOffshrBsline5MW_Onshore_ElastoDyn_8mps.dat\"    EDFile";
       print "\"unused\"      BDBldFile(1)";
       print "\"unused\"      BDBldFile(2)";
       print "\"unused\"      BDBldFile(3)";
       print "\"NRELOffshrBsline5MW_Onshore_ServoDyn_WT1.dat\"    ServoFile";
       ins=1 } { print }' \
    "$FG/farm_multirotor/FFTest_WT1.fst" > "$FG/farm_multirotor/.tmp" \
    && mv "$FG/farm_multirotor/.tmp" "$FG/farm_multirotor/FFTest_WT1.fst"
runfarm farm_multirotor "$NROTOR_MSG" yes

# A mirrored turbine in a farm.  The mirror is confined to its own OpenFAST
# instance and the blade surfaces do render correctly, but the wake coupling has
# never been measured: FWrap_CalcOutput builds the skew angle from a cross
# product of the disk-averaged wind with the disk normal, which is the
# pseudovector pattern that carries a sign at every other module boundary here.
# Until that is checked the combination is refused rather than run quietly.
mkfarm farm_mirror
sed -i -E 's/^([[:space:]]*)F([[:space:]]+MirrorRotor)/\1       True\2/' "$FG/farm_mirror/FFTest_WT1.fst"
runfarm farm_mirror "$MIRROR_MSG" yes
