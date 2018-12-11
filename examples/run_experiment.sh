#!/usr/bin/env bash

# Script to run the experiment for statistics.

N_trials=10

tableMin=0.001
tableMax=30.0
tableTol=1e-2

N_per_trial=1000000
startSeed=19

executable=./experiment

fileName=$1

runLine="${executable} ${tableMin} ${tableMax} ${tableTol} ${N_trials} ${N_per_trial} ${startSeed}"
echo ${runLine}
eval ${runLine}
