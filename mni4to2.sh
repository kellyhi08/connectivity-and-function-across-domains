#!/bin/bash
#converts 4mm data in MNI to 2mm.  specify the paths
#eg mni4to2.sh /path/to/input /path/to/output
# you may need to replace oldtarg & newtarg with the paths to your FSL standard data
in=$1
out=$2

targ=/usr/local/fsl/5.0.10/data/standard/MNI152_T1_2mm.nii.gz # this is the space we want to move the neuroquery maps too
mri_vol2vol --mov $in --targ $targ --regheader --o $out
