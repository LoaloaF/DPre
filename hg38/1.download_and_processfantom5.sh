# 

# Bash script to download and process the FANTOM5 data, into a format suitable for DPre to use.

# First get the FANTOM5 hg38 raw data, if nto alread downloaded:
wget -c http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_ann.txt.gz
wget -c http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz\

# Run the python  script to extract remove the cancerous cell types, add the annotations, 
# merge replicates (where possible), and nicely tidy everything up.

python3 pack_fantom5_data_human.py