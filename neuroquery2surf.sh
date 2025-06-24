#!/bin/bash

in=$1
hcppath=$2
yourpath=$3
wb_command -volume-to-surface-mapping  ${yourpath}/inputs_query_new/${in}_2.nii.gz ${hcppath}/100307/MNINonLinear/fsaverage_LR32k/100307.L.midthickness_MSMAll.32k_fs_LR.surf.gii ${yourpath}/surf_query/$in.L.func.gii -trilinear

wb_command -volume-to-surface-mapping  ${yourpath}/inputs_query_new/${in}_2.nii.gz ${hcppath}/100307/MNINonLinear/fsaverage_LR32k/100307.R.midthickness_MSMAll.32k_fs_LR.surf.gii ${yourpath}/surf_query/$in.R.func.gii -trilinear
