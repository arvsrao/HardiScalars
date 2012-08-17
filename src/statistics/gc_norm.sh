#! /bin/sh

if [ "X$1" != "X" ] ; then
       ID=$1
fi

if [ "X$ID" = "X"  ] ; then
    echo "No ID given" 1>&2
    exit 1
fi

OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}/Norm_FOD

if [[ ! -f ${OutputDir}/${ID}0.nii.gz ]]; then 
      peakFinding -d ${OutputDir}/${ID}_FOD_Warp_Norm.nii.gz -v -p ${ID} -N 0 -n 2 -o ${OutputDir}
fi

if [[ ! -f ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz && ! -f ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz ]]; then 
      RelabelPeaks.py ${ID} 
else 
   echo "One or both relabled peak files exist"      
fi

if [[ ! -f ${OutputDir}/${ID}_relabled_1st_peak_gc.nii.gz && -f ${OutputDir}/${ID}_Norm_FOD_Warp_relabeled_1st_peak.nii.gz ]]; then 
      GeodesicConcentration -a ${RshDir}/${ID}_FOD_Warp_Norm.nii.gz -b ${OutputDir}/${ID}_Norm_FOD_Warp_relabeled_1st_peak.nii.gz -p ${OutputDir}/${ID}_relabeled_1st_peak
fi

if [[ ! -f ${OutputDir}/${ID}_relabed_2nd_peak_gc.nii.gz && -f ${OutputDir}/${ID}_Norm_FOD_Warp_relabeled_2nd_peak.nii.gz ]]; then 
      GeodesicConcentration -a ${RshDir}/${ID}_FOD_Warp_Norm.nii.gz -b ${OutputDir}/${ID}_Norm_FOD_Warp_relabeled_2nd_peak.nii.gz -p ${OutputDir}/${ID}_relabeled_2nd_peak
fi
