#!/usr/bin/env bash
#
# check_tabs.sh -- verify that no Fortran source file contains tab characters.
#
# OpenFAST style requires 3-space indentation; tab characters are not allowed
# because they render inconsistently across editors and break the alignment of
# continuation lines.
#
# Scans every git-tracked *.f90 / *.F90 file.  Registry-generated *_Types.f90
# files are excluded, since they are produced by the OpenFAST Registry rather
# than edited by hand.
#
# Exit status:
#    0  no tab characters found
#    1  one or more tab characters found
#    2  the check could not run (not a git checkout, or a file failed to scan)
#
set -uo pipefail

TAB=$(printf '\t')
MARK='--->'

repo_root=$(git rev-parse --show-toplevel) || {
   echo "ERROR: not inside a git checkout; cannot determine repository root."
   exit 2
}
cd "$repo_root" || exit 2

echo "OpenFAST Fortran source style check"
echo "==================================="
echo
echo "OpenFAST style requires 3 space indentation in all Fortran source files."
echo "Tab characters are not allowed."
echo

mapfile -t FILES < <(git ls-files -- '*.f90' '*.F90' | grep -v '_Types\.f90$' | sort)

if [ "${#FILES[@]}" -eq 0 ]; then
   echo "ERROR: no Fortran source files found -- is this a git checkout?"
   exit 2
fi

n_bad_files=0
n_bad_lines=0

for file in "${FILES[@]}"; do
   # grep exits 0 on a match, 1 on no match, and >1 on an actual error.
   # Only "no match" may be skipped; a scan failure must not look like a pass.
   matches=$(grep -n -- "$TAB" "$file")
   status=$?
   if [ "$status" -eq 1 ]; then
      continue
   elif [ "$status" -ne 0 ]; then
      echo "ERROR: failed to scan ${file} (grep exit ${status})."
      exit 2
   fi

   n_bad_files=$((n_bad_files + 1))
   echo "$file"

   while IFS= read -r match; do
      lineno=${match%%:*}
      content=${match#*:}
      column=$(awk -v s="$content" 'BEGIN { print index(s, "\t") }')
      n_bad_lines=$((n_bad_lines + 1))

      printf '   line %s, column %s: %s\n' "$lineno" "$column" "${content//$TAB/$MARK}"

      # Inline annotation on the pull request diff.
      if [ "${GITHUB_ACTIONS:-}" = "true" ]; then
         echo "::error file=${file},line=${lineno},col=${column}::Tab character found. Use 3 space indentation; tabs are not allowed."
      fi
   done <<< "$matches"

   echo
done

echo "Checked ${#FILES[@]} Fortran source file(s)."

if [ "$n_bad_files" -gt 0 ]; then
   echo
   echo "ERROR: found ${n_bad_lines} line(s) containing tab characters in ${n_bad_files} file(s)."
   echo "       Replace each tab with spaces (3 space indentation) and commit the fix."
   echo "       Tabs above are shown as '${MARK}'."
   exit 1
fi

echo "No tab characters found."
exit 0
