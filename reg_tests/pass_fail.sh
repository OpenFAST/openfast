#!/bin/bash

determine_pass_fail() {
  diffAnywhere=0
  testSolution=$1
  truthSolution=$2
  tolerance=$3

  if [ ! -f $testSolution ]; then
    diffAnywhere=1
  elif [ ! -f $truthSolution ]; then
    diffAnywhere=1
  else
    python `pwd`/../../reg_tests/compareTwoFASTruns.py $testSolution $truthSolution $tolerance
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
  tolerance=$3

  determine_pass_fail $testSolution $truthSolution $tolerance
  passStatus="$?"

  if [ ${passStatus} -eq 0 ]; then
    exit 0
  else
    exit 1
  fi
}

main "$@"
