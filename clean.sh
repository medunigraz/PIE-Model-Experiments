#!/bin/bash

shopt -s extglob

rm -vr !("README.md"|"requirements.txt"|".gitignore"|"clean.sh"|"eval.sh"|".git"|"scripts"|"calibration"|"setups"|"results")

shopt -u extglob
