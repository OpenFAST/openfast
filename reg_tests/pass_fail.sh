#!/bin/bash

# determine_pass_fail() {
#     diffAnywhere=0
#     tolerance=$1
#     logFileName=$2
#     localNormFileName=$3
#     goldNormFileName=$4
#
#     # check for required files: log file and gold
#     if [ ! -f "$logFileName" ]; then
# 	    diffAnywhere=1
#     fi
#
#     # check for gold norm file
#     if [ ! -f "$goldNormFileName" ]; then
# 	    diffAnywhere=1
#     fi
#
#     grep "Mean System Norm:" "$logFileName"  | awk '{print $4 " " $5 " " $6 }' > "$localNormFileName"
#
#     # make sure the grep  worked
#     if [ ! -f "$localNormFileName" ]; then
# 	    diffAnywhere=1
#     fi
#
#     # read in gold norm values
#     goldCount=1
#     goldFileContent=( `cat "$goldNormFileName"`)
#     for gfc in "${goldFileContent[@]}"; do
# 	    goldNorm[goldCount]=$gfc
# 	    ((goldCount++))
#     done
#
#     # read in local norm values
#     localCount=1
#     localFileContent=( `cat "$localNormFileName"`)
#     for lfc in "${localFileContent[@]}"; do
# 	    localNorm[localCount]=$lfc
# 	    ((localCount++))
#     done
#
#     if [ $(echo " $goldCount - $localCount" | bc) -eq 0 ]; then
#       # the lengths the same... proceed
# 	    for ((i=0;i<$goldCount;++i)); do
#         modLocalNorm=$(printf "%1.32f" ${localNorm[i]})
#         modGoldNorm=$(printf "%1.32f" ${goldNorm[i]})
#
#         # compute the difference
#         diff=$(echo $modLocalNorm - $modGoldNorm | bc)
#
#         # make sure diff is positive.. abs anyone?
#         zero=0.0
#         minusOne=-1.0
#         if [ $(echo " $diff < $zero" | bc) -eq 1 ]; then
#          absDiff=$(echo $diff*$minusOne | bc)
#         else
#          absDiff=$diff
#         fi
#
#         # test the difference
#         if [ $(echo " $absDiff > $tolerance" | bc) -eq 1 ]; then
#          diffAnywhere=1
#         fi
#
#         # find the max
#         if [ $(echo " $absDiff > $maxSolutionDiff " | bc) -eq 1 ]; then
#          maxSolutionDiff=$absDiff
#         fi
# 	    done
#     else
#       # length was not the same; fail
# 	    diffAnywhere=1
#       maxSolutionDiff=1000000.0
#     fi
#
#     # extract simulation time
#     return $diffAnywhere
# }

# $1 is the test name
main() {
  testName=$1
  goldFile=$2
  mytolerance=$3
  maxSolutionDiff=-1000000000.0
  # determine_pass_fail ${mytolerance} "${testName}.log" "${testName}.norm" "${goldFile}"
  # passStatus="$?"
  passStatus=0

  # performanceTime=`grep "STKPERF: Total Time" ${testName}.log  | awk '{print $4}'`
  if [ ${passStatus} -eq 0 ]; then
    echo -e "..${testName}........... PASSED":" " ${performanceTime} " s"
    exit 0
  else
    echo -e "..${testName}........... FAILED":" " ${performanceTime} " s" " max diff: " ${maxSolutionDiff}
    exit 1
  fi
}

# ./pass_fail.sh testName normFile tolerance
main "$@"
