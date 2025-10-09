export WIREMOD_WORKING_DIR="/exp/sbnd/app/users/aantonak/WireModProd/SBNDWireModAna"
export CALIB_WORKING_DIR="/exp/sbnd/app/users/aantonak/SBND_calib_recom"
export DATA_PATH=$CALIB_WORKING_DIR/data/
export PLOT_PATH=$CALIB_WORKING_DIR/output/plots/
export OUTPUTROOT_PATH=$CALIB_WORKING_DIR/output/root/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$CALIB_WORKING_DIR/include/:$WIREMOD_WORKING_DIR/include_wire/
#export ROOT_INCLUDE_PATH=$CALIB_WORKING_DIR/include/:$WIREMOD_WORKING_DIR/include/
source $CALIB_WORKING_DIR/bin/BashColorSets.sh

#####################################################################################
## -- Host dependent settings

#### -- Setup root
MY_OS_REL=$(cat /etc/os-release | grep ^NAME | sed -e 's/NAME=//g' -e 's/"//g')
if [[ "$MY_OS_REL" == "AlmaLinux" && $(hostname) != *"dune-gpu01"* ]]; then
  source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
  spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
elif [ "$MY_OS_REL" = "Scientific Linux" ]; then

  #Check if PRODUCTS is undefined -- if so, set up relevant ups area
  if [[ -z $PRODUCTS ]]; then
    if [[ $HOSTNAME == "sbnd"* ]]; then
      echo "SBND"
      MY_EXPERIMENT="sbnd"
    elif [[ $HOSTNAME == "icarusgpvm"* ]]; then
      echo "ICARUS"
      MY_EXPERIMENT="sbnd"
    else
      echo "Warning: Unrecognized hostname $HOSTNAME"
    fi
  fi

  source /cvmfs/${MY_EXPERIMENT}.opensciencegrid.org/products/${MY_EXPERIMENT}/setup_${MY_EXPERIMENT}.sh
  setup root v6_28_12 -q e26:p3915:prof 
  setup xrootd v5_5_5a -q e26:p3915:prof
  setup cmake v3_27_4
else
  echo "WARNING: Seems you are using a private machine to run this repo"
  echo "I do not automatically set up ROOT. If ROOT is already setup, it shoould be okay"
fi

#### -- setup grid output
if [[ $(hostname) == *sbnd* || $(hostname) == *jupyter* ]]; then
    export SBNDCALIB_GRID_OUT_DIR="/pnfs/sbnd/scratch/users/$USER/sbnd_calib_out"
    mkdir -p $SBNDCALIB_GRID_OUT_DIR
fi

#### -- Calib ntuple list dir
export SAMPLE_PATH=$DATA_PATH/sample_list/sungbinosx/
if [[ `hostname` == *"sbnd"* || `hostname` == *"$USER"* ]]; then
    source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
    export SAMPLE_PATH=$DATA_PATH/sample_list/sbndgpvm/
    export SBND_DATA_PATH=/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/
    #export SBNDDATA_VERSION=v01_28_00
    export SBNDDATA_VERSION=v01_30_00

    export SBND_YZCORR_PATH=/exp/sbnd/app/users/yadav/Calibration/sbndcode_v10_06_01/sbnd_data/v01_33_00/YZmaps/
fi

if [[ `hostname` == *"jupyter"* ]]
then
    source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
    export SAMPLE_PATH=$DATA_PATH/sample_list/sbndgpvm/
    export SBND_DATA_PATH=/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/
    #export SBNDDATA_VERSION=v01_28_00
    export SBNDDATA_VERSION=v01_30_00

    export SBND_YZCORR_PATH=/exp/sbnd/data/users/sungbino/sbnd_corr/yzcorr/
    export SBND_YZCORR_VERSION=202506
fi

if [[ `hostname` == *"dune-gpu01"* ]]
then
    export SAMPLE_PATH=$DATA_PATH/sample_list/dune-gpu01/
    export FILELIST_LABEL=_dune_gpu01_run_
fi
#####################################################################################
