#!/bin/bash
# TODO:
# - visualize_calibration.py
# - error_density.py
# - phase_singularitiy_tracking.py
# - ecg_comparison.py
# - ecg_animation.py

# external tools ==============================================================
PIE_EXE=pie-solver
MT_EXE=meshtool
P3_EXE=python3
CARP_EXE=openCARP
GIZMO_EXE=gizmo

# Variables ===================================================================
NP=32
FEM_THRESH="400"
TESTCASES=("A2_functional") # "A1_anatomical" "A2_functional" "A3_wholeheart" "B1_restitution" "B2_curvature" "B3_diffusion"

declare -a RES_RD_UM=("250")
declare -a RES_PIE_UM=("1000")

# Functions ===================================================================
function create_dir {
  if [ ! -d $1 ]; then
    mkdir -p $1
    echo "Created output directory ${1}."
  fi
}

function run_pie {
  local EXP="$1"
  local RES="$2"
  local MSH="$3"
  local PLN="$4"
  local SIM="$5"
  local ENS="$6"
  local PAR="$7"
  local SUF="$8"

  local SIM_DIR_PIE="${SIM}/pie_${RES}um"
  local ENS_DIR_PIE="${ENS}/pie_${RES}um"

  if [ "$SUF" != "" ]; then
    SIM_DIR_PIE="${SIM_DIR_PIE}_${SUF}"
    ENS_DIR_PIE="${ENS_DIR_PIE}_${SUF}"
  fi

  create_dir $SIM_DIR_PIE
  create_dir $ENS_DIR_PIE

  $PIE_EXE --msh=$MSH --pln=$PLN --out="${SIM_DIR_PIE}/" $PAR
  
  # convert time slice data from igb to paraview ens format for visualization
  local NOD_DATA=""

  for ENS_DATA in "${SIM_DIR_PIE}/"*.igb; do
    IGB_FILENAME=$(basename -- $ENS_DATA)
    IGB_FILENAME="${IGB_FILENAME: 0:-4}"

    if [ $IGB_FILENAME == "lat_ms" ] || [ $IGB_FILENAME == "lrt_ms" ]; then
      continue;
    fi

    if [ $IGB_FILENAME == "edx_n" ]; then
      NOD_DATA="${NOD_DATA} -ele=${ENS_DATA}"
    else
      NOD_DATA="${NOD_DATA} -nod=${ENS_DATA}"
    fi
  done
  
  $MT_EXE collect -imsh=$MSH -omsh="${ENS_DIR_PIE}/pie_data" $NOD_DATA -ifmt=carp_bin -ofmt=ens_bin
}

function compare_rd_pie {
  local RES_RD_UM="$1"
  local RES_PIE_UM="$2"
  local SETUP="$3"
  local SIM="$4"
  local ENS="$5"
  local SUF="$6"

  local SIM_COMP="${SIM}/comp_${RES_PIE_UM}"
  local ENS_COMP="${ENS}/comp_${RES_PIE_UM}"
  local MSH_RD="${SETUP}mesh_${RES_RD_UM}um/mesh_${RES_RD_UM}um"
  local MSH_PIE="${SETUP}mesh_${RES_PIE_UM}um/mesh_${RES_PIE_UM}um"
  local SIM_RD="${SIM}/rd_${RES_RD_UM}um/sim"
  local SIM_PIE="${SIM}/pie_${RES_PIE_UM}um"

  if [ "$SUF" != "" ]; then
    SIM_COMP="${SIM_COMP}_${SUF}"
    ENS_COMP="${ENS_COMP}_${SUF}"
    SIM_PIE="${SIM_PIE}_${SUF}"
  fi

  create_dir $SIM_COMP
  create_dir $ENS_COMP

  $PIE_EXE --compare --mshA=$MSH_RD --latA="${SIM_RD}/lats-thresh-2.dat" --lrtA="${SIM_RD}/lrts-thresh-2.dat" --mshB=$MSH_PIE --latB="${SIM_PIE}/lat_all_ms.dat" --lrtB="${SIM_PIE}/lrt_all_ms.dat" --out=$SIM_COMP --nodal --np=$NP
  NOD_DATA=""

  for IGB_FILEPATH in "${SIM_COMP}/"*.igb; do
    IGB_FILENAME=$(basename -- $IGB_FILEPATH)
    IGB_FILENAME="${IGB_FILENAME: 0:-4}"
    NOD_DATA="${NOD_DATA} -nod=${IGB_FILEPATH}"
  done

  $MT_EXE collect -imsh=$MSH_PIE -omsh="${ENS_COMP}/compare.case" $NOD_DATA -ifmt=carp_bin -ofmt=ens_bin
  $P3_EXE ./scripts/error-density.py --compdir=$SIM_COMP --outdir=$SIM
}

# Parse Inputs ====================================================================================
#-exp=...
for i in "$@"
do
case $i in
  -np=*)
  NP="${i#*=}"
  ;;
  *)
  # unknown option
  echo "Unknown parameter passed: $1"; exit 1
  ;;
esac
done

# Calibration =================================================================
if false; then 
  $PIE_EXE --tune --pln=./calibration/forcepss/template.plan.json --out=./calibration/forcepss/tune --verb --np=$NP
  $P3_EXE ./scripts/visualize-calibration.py --pln=./calibration/forcepss/tune/calibrated-functions.plan.json
  #cp -a ./calibration/forcepss/tune/calibration/. ./calibration/
fi

# Evaluate ====================================================================
for TESTCASE in $TESTCASES; do
  SETUP_DIR="./setups/${TESTCASE}/"
  SIM_DIR="./results/sim/${TESTCASE}"
  ENS_DIR="./results/ens/${TESTCASE}"

  for PLN_ID in $SETUP_DIR*.json; do
    PLN_PATH="${PLN_ID}" || break
  done
  
  create_dir $SIM_DIR
  create_dir $ENS_DIR

  # Launch RD Simulations -----------------------------------------------------
  if true; then
    for RES_UM in "${RES_RD_UM[@]}"; do
      if [ $RES_UM -gt $FEM_THRESH ]; then
        continue;
      fi

      MSH_DIR="${SETUP_DIR}/mesh_${RES_UM}um"
      MSH_ID=$(basename -- $MSH_DIR)
      MSH_PATH="${MSH_DIR}/${MSH_ID}"
      SIM_RD="${SIM_DIR}/rd_${RES_UM}um"
      ENS_RD="${ENS_DIR}/rd_${RES_UM}um"

      if [ ! -d $MSH_DIR ]; then
        $P3_EXE ./scripts/meshgen.py --msh=$TESTCASE --res=$RES_UM --out=$SETUP_DIR
        $MT_EXE convert -imsh=$MSH_PATH -ifmt=carp_bin -omsh=$MSH_PATH -ofmt vtk_bin
      fi

      #if [ ! -d $SIM_RD ]; then
        create_dir $SIM_RD

        $PIE_EXE --pln2par --msh=$MSH_PATH --pln=$PLN_PATH --out=$SIM_RD --verb
        mpirun -np $NP $CARP_EXE +F "${SIM_RD}/lat.par"
        mpirun -np $NP $CARP_EXE +F "${SIM_RD}/prp.par"
        mpirun -np $NP $CARP_EXE +F "${SIM_RD}/sim.par"

        $P3_EXE ./scripts/apply_tstart_offset.py --msh=$MSH_PATH --dat="${SIM_RD}/sim/lats-thresh.dat" --offset=1800 --out="${SIM_RD}/sim/lats-thresh-2.dat"
        $P3_EXE ./scripts/apply_tstart_offset.py --msh=$MSH_PATH --dat="${SIM_RD}/sim/lrts-thresh.dat" --offset=1800 --out="${SIM_RD}/sim/lrts-thresh-2.dat"
      #fi

      #if [ ! -d $ENS_RD ]; then
        create_dir $ENS_RD
        $MT_EXE collect -imsh=$MSH_PATH -omsh="${ENS_RD}/rd_data" -nod="${SIM_RD}/sim/vm.igb" -ifmt=carp_bin -ofmt=ens_bin
      #fi
    done
  fi

  # Launch PIE Model Simulations ------------------------------------------------------------------
  if true; then 
    for RES_UM in "${RES_PIE_UM[@]}"; do
      MSH_DIR="${SETUP_DIR}/mesh_${RES_UM}um"
      MSH_ID=$(basename -- $MSH_DIR)
      MSH_PATH="${MSH_DIR}/${MSH_ID}"

      if [ ! -d $MSH_DIR ]; then
        $P3_EXE ./scripts/meshgen.py --msh=$TESTCASE --res=$RES_UM --out=$SETUP_DIR
        $MT_EXE convert -imsh=$MSH_PATH -ifmt=carp_bin -omsh=$MSH_PATH -ofmt vtk_bin
      fi
      
      if [ "$TESTCASE" = "A1_anatomical" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --dat=32"
      elif [ "$TESTCASE" = "A2_functional" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --dat=32"
      elif [ "$TESTCASE" = "A3_wholeheart" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --dat=32"
        $GIZMO_EXE --model-name $MSH_PATH --model-format carp_bin --sim-plan $PLN_PATH --extern-vm "/home/tom/workspace/eikonal-experiments/results/sim/7_whole-heart/pie_2000um/vm_mv.igb" --trace-output json --lin-solver amgcl --output-dir 'wholeheartecgs' --verbose --np $NP
      elif [ "$TESTCASE" = "B1_restitution" ]; then
        #run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-cvr --no-apdr --no-curv --no-dif --dat=32" "00"
        #run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-cvr --no-curv --no-dif --dat=32" "01"
        #run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-apdr --no-curv --no-dif --dat=32" "10"
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-curv --no-dif --dat=32 --prep" "11"
      elif [ "$TESTCASE" = "B2_curvature" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-dif --dat=32"
      elif [ "$TESTCASE" = "B3_diffusion" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-curv --dat=32"
      else
        echo "Error: invalid testcase specified, must be one of the following:"
        echo "(A1_anatomical, A2_functional, A3_wholeheart, B1_restitution, B2_curvature, B3_diffusion)"
      fi
    done
  fi

  # Compare RD and PIE Results -------------------------------------------------
  if false; then
    for RES_RD in "${RES_RD_UM[@]}"; do
      for RES_PIE in "${RES_PIE_UM[@]}"; do
        if   [ "$TESTCASE" = "B1_restitution" ]; then
          compare_rd_pie $RES_RD $RES_PIE $SETUP_DIR $SIM_DIR $ENS_DIR "11"
        elif [ "$TESTCASE" = "B2_curvature" ] || [ "$TESTCASE" = "B3_diffusion" ]; then
          compare_rd_pie $RES_RD $RES_PIE $SETUP_DIR $SIM_DIR $ENS_DIR
        elif [ "$TESTCASE" = "A2_functional" ]; then
          $P3_EXE ./scripts/phase_singularitiy_tracking.py
        elif [ "$TESTCASE" = "A3_wholeheart" ]; then
          $P3_EXE ./scripts/ecg_comparison.py
          $P3_EXE ./scripts/ecg_animation.py
        fi
      done
    done
  fi
done

