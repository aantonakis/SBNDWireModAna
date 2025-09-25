# Setup grid submission

echo "@@hostname: ${hostname}"
hostname

outDir=$1
echo "@@ outDir : ${outDir}"
DFPREFIX=$2

nProcess=$PROCESS
echo "@@ nProcess : "${nProcess}

echo "@@ pwd"
pwd
echo "@@ ls -alh"
ls -alh

thisOutputCreationDir=`pwd`
filesFromSender=${CONDOR_DIR_INPUT}/bin_dir/

echo "@@ check filesFromSender dir"
ls -alh ${filesFromSender}

echo "@@ copy scripts"
mkdir -p data/sample_list/sbndgpvm/
cp ${filesFromSender}/input_list_${nProcess}.txt ./data/sample_list/sbndgpvm/
ls -alh ./data/sample_list/sbndgpvm/
cp ${filesFromSender}/setup_grid.sh .
cp -r ${filesFromSender}/bin .
cp -r ${filesFromSender}/include .
cp -r ${filesFromSender}/include_wire .
cp -r ${filesFromSender}/ndhist_charges_tpc_crossers_grid_ntrack.C .

echo "@@ make output/root"
mkdir -p output/root
#ifdh  mkdir_p ${outDir}
ls -alh

echo "@@ source setup_grid.sh"
source setup_grid.sh
export SBND_YZCORR_PATH=${filesFromSender}/
spack find ifdhc
#spack load ifdhc@2.8.0 arch=linux-almalinux9-x86_64_v3
#spack load ifdhc@2.7.2
spack load ifdhc@2.6.20

spack find root

echo "@@ run"
# Booleans: apply_sce, apply_yz, apply_elife, apply_recomb, IsData

# //////////////////////////////////////////////////////////


# All calibrations for DATA
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, true)" &> log_${nProcess}.log

# SCE Only
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, false, false, false, true)" &> log_${nProcess}.log

# YZ only
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", false, true, false, false, true)" &> log_${nProcess}.log

# E-life only 
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", false, false, true, false, true)" &> log_${nProcess}.log


# ////////////// Simulation Jobs ///////////////////////////

# all calibrations for mc
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, false)" &> log_${nProcess}.log

# No calibrations
root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", false, false, false, false, false)" &> log_${nProcess}.log

# SCE Only
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, false, false, false, false)" &> log_${nProcess}.log

# YZ only
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", false, true, false, false, false)" &> log_${nProcess}.log

# E-life only 
#root -l -b -q "ndhist_charges_tpc_crossers_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", false, false, true, false, false)" &> log_${nProcess}.log


# //////////////////////////////////////////////////////////

ls -alh

echo "@@ setup IFDH_CP_MAXRETRIES to 2"
export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir}

outFILE=${thisOutputCreationDir}/output/root/output_ndhist_charges_tpc_crossers_${nProcess}.root
#outFILE=${thisOutputCreationDir}/${outDir}/output_ndhist_charges_tpc_crossers_${nProcess}.root
if [ -f "$outFILE" ]; then
  echo "ifdh cp ${thisOutputCreationDir}/output/root/output_ndhist_charges_tpc_crossers_${nProcess}.root ${outDir}/${DFPREFIX}_${nProcess}.root"
  ifdh cp ${thisOutputCreationDir}/output/root/output_ndhist_charges_tpc_crossers_${nProcess}.root ${outDir}/${DFPREFIX}_${nProcess}.root
  echo "ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log"
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "@@ Done!"
else
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "File not exist"
fi

