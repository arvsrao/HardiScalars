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

if [[ -f ${OutputDir}/${ID}0.nii.gz && ${OutputDir}/${ID}1.nii.gz ]]; then 
      mv ${OutputDir}/${ID}0.nii.gz ${OutputDir}/${ID}_Norm_FOD_Warp_first_peak.nii.gz
      mv ${OutputDir}/${ID}1.nii.gz ${OutputDir}/${ID}_Norm_FOD_Warp_second_peak.nii.gz
else
   echo "Peakfinding Failed!!"
   exit 1
fi

if [[ ! -f ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz && ! -f ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz ]]; then 
      RelabelPeaks.py ${ID} 
else 
   echo "One or both relabled peak files exist!!"      
fi

if [[ ! -f ${OutputDir}/${ID}_relabeled_1st_peak_gc.nii.gz && -f ${OutputDir}/${ID}_FOD_Warp_Norm.nii.gz ]]; then 
      GeodesicConcentration -a ${OutputDir}/${ID}_FOD_Warp_Norm.nii.gz -b ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz -p ${OutputDir}/${ID}_relabeled_1st_peak
fi

if [[ ! -f ${OutputDir}/${ID}_relabeled_2nd_peak_gc.nii.gz && -f ${OutputDir}/${ID}_FOD_Warp_Norm.nii.gz ]]; then 
      GeodesicConcentration -a ${OutputDir}/${ID}_FOD_Warp_Norm.nii.gz -b ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz -p ${OutputDir}/${ID}_relabeled_2nd_peak
fi

#compute peak magnitudes.
if [[ ! -f ${OutputDir}/${ID}_relabeled_1st_pm.nii.gz ]]; then 
      3dcalc -a ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz[0] -b ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz[1] -c ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz[2] -expr 'sqrt(a*a + b*b + c*c)' -session ${OutputDir} -prefix ${ID}_relabeled_1st_pm.nii.gz
fi

if [[ ! -f ${OutputDir}/${ID}_relabeled_2nd_pm.nii.gz ]]; then 
      3dcalc -a ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz[0] -b ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz[1] -c ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz[2] -expr 'sqrt(a*a + b*b + c*c)' -session ${OutputDir} -prefix ${ID}_relabeled_2nd_pm.nii.gz
fi

#compute sum measure
if [[ ! -f ${OutputDir}/${ID}_sum_measure.nii.gz ]]; then 
      3dcalc -a ${OutputDir}/${ID}_relabeled_1st_pm.nii.gz -b ${OutputDir}/${ID}_relabeled_1st_peak_gc.nii.gz -c ${OutputDir}/${ID}_relabeled_2nd_peak_gc.nii.gz -expr 'a*b + (1-a)*c' -session ${OutputDir} -prefix ${ID}_sum_measure.nii.gz 
fi
 
