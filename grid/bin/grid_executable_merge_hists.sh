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
cp -r ${filesFromSender}/merge_hists_grid.C .

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


# inputs: apply_sce, apply_yz, apply_elife, apply_recomb, IsData, dim, tracks_sel, crt_sel, pathological_sel, lifetime_sel

# dim = {x, y, z, txz, tyz, dq/dx, Q, width, goodness, pathological}

# Dimensions
# 0 = x
# 1 = y
# 2 = z
# 3 = txw
# 4 = tyz
# 5 = dQdx
# 6 = Q
# 7 = width
# 8 = Goodness
# 9 = Pathological

root -l -b -q "merge_hists_grid.C(\"input_list_${nProcess}.txt\", \"${nProcess}\", {0, 1, 2})" &> log_${nProcess}.log


# //////////////////////////////////////////////////////////

ls -alh

echo "@@ setup IFDH_CP_MAXRETRIES to 2"
export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir}

outFILE=${thisOutputCreationDir}/output/root/output_merge_hists_${nProcess}.root
#outFILE=${thisOutputCreationDir}/${outDir}/output_ndhist_charges_tracks_crossers_${nProcess}.root
if [ -f "$outFILE" ]; then
  echo "ifdh cp ${thisOutputCreationDir}/output/root/output_merge_hists_${nProcess}.root ${outDir}/${DFPREFIX}_${nProcess}.root"
  ifdh cp ${thisOutputCreationDir}/output/root/output_merge_hists_${nProcess}.root ${outDir}/${DFPREFIX}_${nProcess}.root
  echo "ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log"
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "@@ Done!"
else
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "File not exist"
fi

