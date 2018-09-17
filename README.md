This repository is a working repository for physics analysis in the CMS experiment, EXO group, MET+X subgroup, Semi-Visible Jets analysis.

Most of what is here is to be used with CMSSW_8_0_28. No guarantee for support is made.


Current usage:
run in the directory with main.py
macro_name.py must be in /macros/
ROOT_file_list.txt must be in /input_conf/
TTree_name_list.txt must be in /input_conf/
localStorageDirectory must exist
outputFileName.root will be called with 'RECREATE'

python main.py [macro_name.py] [ROOT_file_list.txt] [TTree_name_list.txt] [localStorageDirectory] [outputFileName.root]

reccomended to run macros/makeFriend.py first. This macro create a 'friend tree' with information regarding assigning GenParticles to certain AK8Jets.
