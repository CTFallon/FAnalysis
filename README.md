This repository is a working repository for physics analysis in the CMS experiment, EXO group, MET+X subgroup, Semi-Visible Jets analysis.

Most of what is here is to be used with CMSSW_8_0_28. No garuntee for support (nor spelling) is made.

Many things are not optimized (hard-coded file paths, hard-coded parameter values, etc.) and there is currectly no seperation between the code that computes labor-intensive values (i.e., which particles are in jets, how many particles each jet has, whether to include a particle in the jet due to pT values, which jets to include in MT calcualtion....) and the code that makes the plots. Thus, a small change in the plots requires re-computation of all that stuff. First step for improvement is to seperate these bits of code, probably by using the TTree friend mechanism in ROOT to avoid saving the nTuples again....


Current usage:
run in the directory with main.py
macro_name.py must be in /macros/
ROOT_file_list.txt must be in /input_conf/
TTree_name_list.txt must be in /input_conf/
localStorageDirectory must exist
outputFileName.root will be called with 'RECREATE'

pyton main.py <macro_name.py> <ROOT_file_list.txt> <TTree_name_list.txt> <localStorageDirectory> <outputFileName.root>


