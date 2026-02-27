#!/bin/bash

# Executables =====================================================================================
PIE_EXE=./bin/pie-solver
MT_EXE=meshtool
P3_EXE=python3
CARP_EXE=openCARP.opt
GIZMO_EXE=./bin/gizmo

# Default Variables ===============================================================================
NP=32
N_RUNS=5

declare -a TESTCASES=("B1_restitution" "B2_curvature" "B3_diffusion" "A1_anatomical" "A2_functional" "A3_wholeheart")
declare -a RES_RD_UM=("250")
declare -a RES_PIE_UM=("1000")

# Functions =======================================================================================
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
  $P3_EXE ./scripts/error_density.py --compdir=$SIM_COMP --outdir=$SIM
}

# Parse Inputs ====================================================================================
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

# Calibration =====================================================================================
if true; then 
  $PIE_EXE --tune --pln=./calibration/forcepss/template.plan.json --out=./calibration/forcepss/tune --np=$NP
  $P3_EXE ./scripts/visualize_calibration.py --pln=./calibration/forcepss/tune/calibrated-functions.plan.json --odir=./calibration/
fi

# Evaluate ========================================================================================
for TESTCASE in "${TESTCASES[@]}"; do
  SETUP_DIR="./setups/${TESTCASE}/"
  SIM_DIR="./results/sim/${TESTCASE}"
  ENS_DIR="./results/ens/${TESTCASE}"
  PNG_DIR="./results/png/${TESTCASE}"

  for PLN_ID in $SETUP_DIR*.json; do
    PLN_PATH="${PLN_ID}" || break
  done
  
  create_dir $SIM_DIR
  create_dir $ENS_DIR
  create_dir $PNG_DIR

  # Launch RD Simulations -------------------------------------------------------------------------
  if true; then
    for RES_UM in "${RES_RD_UM[@]}"; do
      MSH_DIR="${SETUP_DIR}/mesh_${RES_UM}um"
      MSH_ID=$(basename -- $MSH_DIR)
      MSH_PATH="${MSH_DIR}/${MSH_ID}"
      SIM_RD="${SIM_DIR}/rd_${RES_UM}um"
      ENS_RD="${ENS_DIR}/rd_${RES_UM}um"

      if [ ! -d $MSH_DIR ]; then
        $P3_EXE ./scripts/meshgen.py --msh=$TESTCASE --res=$RES_UM --out=$SETUP_DIR
        $MT_EXE convert -imsh=$MSH_PATH -ifmt=carp_bin -omsh=$MSH_PATH -ofmt vtk_bin
      fi

      if [ ! -d $SIM_RD ]; then
        create_dir $SIM_RD

        if [ "$TESTCASE" = "A3_wholeheart" ]; then
          $PIE_EXE --pln2par --msh=$MSH_PATH --pln="./setups/A3_wholeheart/A3_wholeheart-rd.plan.json" --out=$SIM_RD
        else
          $PIE_EXE --pln2par --msh=$MSH_PATH --pln=$PLN_PATH --out=$SIM_RD
        fi

        mpirun -np $NP $CARP_EXE +F "${SIM_RD}/lat.par"
        mpirun -np $NP $CARP_EXE +F "${SIM_RD}/prp.par"
        mpirun -np $NP $CARP_EXE +F "${SIM_RD}/sim.par"

        $P3_EXE ./scripts/apply_tstart_offset.py --msh=$MSH_PATH --dat="${SIM_RD}/sim/lats-thresh.dat" --offset=1800 --out="${SIM_RD}/sim/lats-thresh-2.dat"
        $P3_EXE ./scripts/apply_tstart_offset.py --msh=$MSH_PATH --dat="${SIM_RD}/sim/lrts-thresh.dat" --offset=1800 --out="${SIM_RD}/sim/lrts-thresh-2.dat"
      fi

      if [ ! -d $ENS_RD ]; then
        create_dir $ENS_RD
        $MT_EXE collect -imsh=$MSH_PATH -omsh="${ENS_RD}/rd_data" -nod="${SIM_RD}/sim/vm.igb" -ifmt=carp_bin -ofmt=ens_bin
      fi
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

      # output map options
      # 1: <reserved>
      # 2: tAp [ms]
      # 3: <reserved>
      # 4: APD_90 [ms]
      # 5: DI [ms]
      # 6: Vm [mV]
      # 7: FSM states [0-8]
      # 8: ta [ms]

      # output options are encoded as 8-bit:
      #      8765 4321
      # 0   (0000 0000): no output
      # 32  (0010 0000): Vm only (default)
      # 96  (0110 0000): Vm and FSM states
      # 224 (1110 0000): Vm, FSM states and ta
      # 255 (1111 1111): all

      # LAT and LRT maps are always generated

      if [ "$TESTCASE" = "A1_anatomical" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --dat=32"
      elif [ "$TESTCASE" = "A2_functional" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --dat=32"
      elif [ "$TESTCASE" = "A3_wholeheart" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --dat=32"
      elif [ "$TESTCASE" = "B1_restitution" ]; then
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-cvr --no-apdr --no-curv --no-dif --dat=32" "00"
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-cvr --no-curv --no-dif --dat=32" "01"
        run_pie $TESTCASE $RES_UM $MSH_PATH $PLN_PATH $SIM_DIR $ENS_DIR "--verb --np=$NP --no-apdr --no-curv --no-dif --dat=32" "10"
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

  # Compare RD and PIE Results --------------------------------------------------------------------
  if true; then
    for RES_RD in "${RES_RD_UM[@]}"; do
      for RES_PIE in "${RES_PIE_UM[@]}"; do
        if   [ "$TESTCASE" = "B1_restitution" ]; then
          compare_rd_pie $RES_RD $RES_PIE $SETUP_DIR $SIM_DIR $ENS_DIR "11"
        elif [ "$TESTCASE" = "B2_curvature" ] || [ "$TESTCASE" = "B3_diffusion" ]; then
          compare_rd_pie $RES_RD $RES_PIE $SETUP_DIR $SIM_DIR $ENS_DIR
        elif [ "$TESTCASE" = "A2_functional" ]; then
          $P3_EXE ./scripts/track_phase_singularity.py --idir=$SIM_DIR --model=rd  --odir=./results/png/A2_functional
          $P3_EXE ./scripts/track_phase_singularity.py --idir=$SIM_DIR --model=pie --odir=./results/png/A2_functional
          $P3_EXE ./scripts/phase_singularity_error.py --simdir=$SIM_DIR --outdir=$SIM_DIR
          #compare_rd_pie $RES_RD $RES_PIE $SETUP_DIR $SIM_DIR $ENS_DIR
          $P3_EXE ./scripts/frequency_maps.py --simdir=$SIM_DIR --outdir=$SIM_DIR
          $P3_EXE ./scripts/phasefield_correlation.py --simdir=$SIM_DIR --outdir=$SIM_DIR
          #PHASE_DIR="${ENS_DIR}/phasemaps"
          #create_dir $PHASE_DIR
          #$MT_EXE collect -imsh="${SETUP_DIR}/mesh_250um/mesh_250um" -omsh="${PHASE_DIR}/phasemap_rd" -nod="${SIM_DIR}/rd_250um/sim/phasemap.igb" -ifmt=carp_bin -ofmt=ens_bin
          #$MT_EXE collect -imsh="${SETUP_DIR}/mesh_1000um/mesh_1000um" -omsh="${PHASE_DIR}/phasemap_pie" -nod="${SIM_DIR}/pie_1000um/phasemap.igb" -ifmt=carp_bin -ofmt=ens_bin
        elif [ "$TESTCASE" = "A3_wholeheart" ]; then
          #$P3_EXE ./scripts/ecg_comparison.py ./results/ecg/ecg_rd_250um.json ./results/ecg/ecg_pie_1000um.json
          $P3_EXE ./scripts/ecg_comparison_v2.py ./results/ecg/ecg_rd_250um.json ./results/ecg/ecg_pie_1000um.json
          $P3_EXE ./scripts/ecg_quantitative.py --ecgA=./results/ecg/ecg_rd_250um.json --ecgB=./results/ecg/ecg_pie_1000um.json --out=./results/A3_wholeheart
        fi
      done
    done
  fi
done

# Outliers and Performance Scaling ================================================================
if true; then
  $P3_EXE ./scripts/plot_outliers.py --outdir=./results

  PERF_DIR="./results/sim/A2c_scaling"
  create_dir $PERF_DIR

  $P3_EXE ./scripts/performance_scaling.py --mode RD,PIE --res 250 --dstart 250 --dend 7000 --dstep 250 --np $NP --plot --outdir $PERF_DIR
fi

# Supplementary Experiments =======================================================================

# A2_variability - Rotor Sensitivity
if true; then
  $PIE_EXE --msh="./setups/A2_functional/mesh_1000um/mesh_1000um" --pln="./setups/A2_functional/A2a_functional.plan.json" --out="./results/sim/A2a_functional/pie_1000um/" --verb --np=$NP --dat=32
  $PIE_EXE --msh="./setups/A2_functional/mesh_1000um/mesh_1000um" --pln="./setups/A2_functional/A2b_functional.plan.json" --out="./results/sim/A2b_functional/pie_1000um/" --verb --np=$NP --dat=32

  $MT_EXE collect -imsh="./setups/A2_functional/mesh_1000um/mesh_1000um" -omsh="./results/ens/A2a_functional/pie_1000um/pie_data" -nod="./results/sim/A2a_functional/pie_1000um/vm_mv.igb" -ifmt=carp_bin -ofmt=ens_bin
  $MT_EXE collect -imsh="./setups/A2_functional/mesh_1000um/mesh_1000um" -omsh="./results/ens/A2b_functional/pie_1000um/pie_data" -nod="./results/sim/A2b_functional/pie_1000um/vm_mv.igb" -ifmt=carp_bin -ofmt=ens_bin
fi

# A3_wholeheart - Impact of Spatial Resampling
if true; then
  # non-preserved
  $P3_EXE "./scripts/performance_eval.py" --msh="./setups/A3_wholeheart/meshes_non_preserved/mesh_1000um/mesh_1000um" --pln="./setups/A3_wholeheart/meshes_non_preserved/A3_wholeheart-1k.plan.json" --out="./results/sim/A3_wholeheart/A1" --np=$NP --runs=$N_RUNS
  $P3_EXE "./scripts/performance_eval.py" --msh="./setups/A3_wholeheart/meshes_non_preserved/mesh_2000um/mesh_2000um" --pln="./setups/A3_wholeheart/meshes_non_preserved/A3_wholeheart-2k.plan.json" --out="./results/sim/A3_wholeheart/A2" --np=$NP --runs=$N_RUNS
  $P3_EXE "./scripts/performance_eval.py" --msh="./setups/A3_wholeheart/meshes_non_preserved/mesh_3000um/mesh_3000um" --pln="./setups/A3_wholeheart/meshes_non_preserved/A3_wholeheart-3k.plan.json" --out="./results/sim/A3_wholeheart/A3" --np=$NP --runs=$N_RUNS

  $P3_EXE "./scripts/ecg_quantitative.py" --ecgA="./results/ecg/ecg_rd_250um.json" --ecgB="./results/sim/A3_wholeheart/A1/ecg.json" --out="./results/A3_wholeheart"
  $P3_EXE "./scripts/ecg_quantitative.py" --ecgA="./results/ecg/ecg_rd_250um.json" --ecgB="./results/sim/A3_wholeheart/A2/ecg.json" --out="./results/A3_wholeheart"
  $P3_EXE "./scripts/ecg_quantitative.py" --ecgA="./results/ecg/ecg_rd_250um.json" --ecgB="./results/sim/A3_wholeheart/A3/ecg.json" --out="./results/A3_wholeheart" # no sustained reentry!

  # surf-preserved
  $P3_EXE "./scripts/performance_eval.py" --msh="./setups/A3_wholeheart/meshes_surf_preserved/mesh_1000um/mesh_1000um" --pln="./setups/A3_wholeheart/meshes_surf_preserved/A3_wholeheart-1k.plan.json" --out="./results/sim/A3_wholeheart/B1" --np=$NP --runs=$N_RUNS
  $P3_EXE "./scripts/performance_eval.py" --msh="./setups/A3_wholeheart/meshes_surf_preserved/mesh_2000um/mesh_2000um" --pln="./setups/A3_wholeheart/meshes_surf_preserved/A3_wholeheart-2k.plan.json" --out="./results/sim/A3_wholeheart/B2" --np=$NP --runs=$N_RUNS
  $P3_EXE "./scripts/performance_eval.py" --msh="./setups/A3_wholeheart/meshes_surf_preserved/mesh_3000um/mesh_3000um" --pln="./setups/A3_wholeheart/meshes_surf_preserved/A3_wholeheart-3k.plan.json" --out="./results/sim/A3_wholeheart/B3" --np=$NP --runs=$N_RUNS

  $P3_EXE "./scripts/ecg_quantitative.py" --ecgA="./results/ecg/ecg_rd_250um.json" --ecgB="./results/sim/A3_wholeheart/B1/ecg.json" --out="./results/A3_wholeheart"
  $P3_EXE "./scripts/ecg_quantitative.py" --ecgA="./results/ecg/ecg_rd_250um.json" --ecgB="./results/sim/A3_wholeheart/B2/ecg.json" --out="./results/A3_wholeheart"
  $P3_EXE "./scripts/ecg_quantitative.py" --ecgA="./results/ecg/ecg_rd_250um.json" --ecgB="./results/sim/A3_wholeheart/B3/ecg.json" --out="./results/A3_wholeheart"

  # visualization
  $P3_EXE ./scripts/ecg_comparison_v2.py ./results/ecg/ecg_rd_250um.json ./results/ecg/ecg_pie_1000um.json ./results/ecg/ecg_pie_1000um_np.json ./results/ecg/ecg_pie_2000um_np.json ./results/ecg/ecg_pie_3000um_np.json ./results/ecg/ecg_pie_1000um_sp.json ./results/ecg/ecg_pie_2000um_sp.json ./results/ecg/ecg_pie_3000um_sp.json
fi

