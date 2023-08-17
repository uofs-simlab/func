#!/usr/bin/env bash

# Script to run experiment.cpp
impl="NonUniformLinearInterpTable<1>"
tableMin=0.001
tableMax=30.0
step=0.1
executable=./plot_impl

runLine="${executable} ${impl} ${tableMin} ${tableMax} ${step}"
echo ${runLine}
eval ${runLine}
