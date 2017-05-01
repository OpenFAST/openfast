#!/bin/bash

determine_pass_fail() {
  diffAnywhere=0
  logFileName=$1
  goldLogFileName=$2

  if [ ! -f ${logFileName} ]; then
    diffAnywhere=1
  elif [ ! -f ${goldLogFileName} ]; then
    diffAnywhere=1
  else
    python `pwd`/../../reg_tests/compareTwoFASTruns.py $logFileName $goldLogFileName
    if [ $? -ne 0 ]; then
      diffAnywhere=1
    fi
  fi
  return $diffAnywhere
}

main() {
  testName=$1
  testSolution=$2
  truthSolution=`pwd`/../../reg_tests/test_files/$testName/$testName.outb

  determine_pass_fail $testSolution $truthSolution
  passStatus="$?"

  if [ ${passStatus} -eq 0 ]; then
    echo -e "..${testName}........... PASSED":" " ${performanceTime} " s"
    exit 0
  else
    echo -e "..${testName}........... FAILED":" " ${performanceTime} " s" " max diff: " ${maxSolutionDiff}
    exit 1
  fi
}

main "$@"
