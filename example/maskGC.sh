#! /bin/bash

# Compute the mean of scalar maps (in this case GFA)
# in a white matter mask. 
#
# 


base_dir=/sbia/home/raoarvin/comp_space/autism_maps
mask=/sbia/home/raoarvin/comp_space/autism_maps/R0090/Norm_FOD/R0090_wmMask.nii.gz

for id in ${base_dir}/* ; do 
      dir=${id##*/}

      if [[ -d ${base_dir}/${dir} ]]; then 
          input=${base_dir}/${dir}/Norm_FOD/${dir}_gfa.nii.gz
          printf "$dir "
          3dBrickStat -mask $mask -mean $input
      fi
done
