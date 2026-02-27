#!/bin/bash

shopt -s extglob

TUNING_DIR="./calibration/forcepss/tune"
SIM_DIR="./results/sim"
ENS_DIR="./results/ens"
ENS_DIR="./results/png"

if [ -d $TUNING_DIR ]; then
  rm -vr $TUNING_DIR
  echo "Removed tuning directory ..."
fi

if [ -d $SIM_DIR ]; then
  rm -vr $SIM_DIR
  echo "Removed simulation directory ..."
fi

if [ -d $ENS_DIR ]; then
  rm -vr $ENS_DIR
  echo "Removed ensight directory ..."
fi

if [ -d $PNG_DIR ]; then
  rm -vr $PNG_DIR
  echo "Removed image directory ..."
fi

#rm -vr !("README.md"|"requirements.txt"|".gitignore"|"clean.sh"|"eval.sh"|".git"|"scripts"|"calibration"|"setups"|"results")

shopt -u extglob
