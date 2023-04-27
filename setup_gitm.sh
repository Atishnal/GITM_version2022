#!/bin/bash
git clone https://github.com/gwbowden/GITM.git
cd GITM
git checkout remotes/origin/unsw_space_chem_all_rr_revised_check
git switch -c remotes/origin/unsw_space_chem_all_rr_revised_check
git pull origin unsw_space_chem_all_rr_revised_check
module load use.own
module load GITM
module load openmpi
module load netcdf
module load openmpi hdf5 intel-compiler

./Config.pl -install -compiler=ifortmpif90 -earth

make          