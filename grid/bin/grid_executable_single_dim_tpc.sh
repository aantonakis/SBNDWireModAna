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
cp -r ${filesFromSender}/single_dim_tpc_grid.C .

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
# inputs: apply_sce, apply_yz, apply_elife, apply_recomb, IsData, dim

# ////////////// Simulation Jobs ///////////////////////////

root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", false, false, false, false, false, 0)" &> log_${nProcess}.log

# all calibrations for mc and x
#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, false, 0)" &> log_${nProcess}.log

#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, false, 1)" &> log_${nProcess}.log

#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, false, 2)" &> log_${nProcess}.log

#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, false, 3)" &> log_${nProcess}.log

#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, false, 4)" &> log_${nProcess}.log


# DATA
#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, true, 0)" &> log_${nProcess}.log

#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, true, 1)" &> log_${nProcess}.log

#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, true, 2)" &> log_${nProcess}.log

#root -l -b -q "single_dim_tpc_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", true, true, true, true, true, 3)" &> log_${nProcess}.log

# //////////////////////////////////////////////////////////

ls -alh

echo "@@ setup IFDH_CP_MAXRETRIES to 2"
export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir}

outFILE=${thisOutputCreationDir}/output/root/output_single_dim_tpc_${nProcess}.root
#outFILE=${thisOutputCreationDir}/${outDir}/output_ndhist_charges_tpc_crossers_${nProcess}.root
if [ -f "$outFILE" ]; then
  echo "ifdh cp ${thisOutputCreationDir}/output/root/output_single_dim_tpc_${nProcess}.root ${outDir}/${DFPREFIX}_${nProcess}.root"
  ifdh cp ${thisOutputCreationDir}/output/root/output_single_dim_tpc_${nProcess}.root ${outDir}/${DFPREFIX}_${nProcess}.root
  echo "ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log"
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "@@ Done!"
else
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "File not exist"
fi

