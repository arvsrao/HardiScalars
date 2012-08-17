#! /bin/sh

if [ "X$1" != "X" ] ; then
       ID=$1
fi

if [ "X$ID" = "X"  ] ; then
    echo "No ID given" 1>&2
    exit 1
fi

OutputDir=/sbia/home/raoarvin/comp_space/autism_maps/${ID}
TplDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/R0090/reg2R0090/fodReg
TplOutputDir=/sbia/home/raoarvin/comp_space/autism_maps/R0090/Norm_FOD

#OutputDir=/sbia/home/raoarvin/comp_space/autism_maps/${ID}/FOD
#TplDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/R0090/reg2R0090/fodReg
#TplOutputDir=/sbia/home/raoarvin/comp_space/autism_maps/R0090/FOD

# Choose the correct path for the Subject based on ID's
if [[ -d  /sbia/sbiaprj/autism/TimRoberts/Processed_HARDI_RAS/${ID}/reg2R0090/ ]]; then  #Tim's dataset
      InputDir=/sbia/sbiaprj/autism/TimRoberts/Processed_HARDI_RAS/${ID}/reg2R0090/fodReg
else 
      InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}/reg2R0090/fodReg      
fi

if [[ ! -d ${OutputDir} ]]; then 
    echo "creating output directory for subject " ${ID}
    mkdir ${OutputDir}
fi

#no matter what create the output directory
OutputDir=/sbia/home/raoarvin/comp_space/autism_maps/${ID}/Norm_FOD
if [[ ! -d ${OutputDir} ]]; then 
    echo "creating output directory for subject " ${ID}
    mkdir ${OutputDir}
fi

if [[ ! -f ${TplOutputDir}/R00900.nii.gz ]]; then
      peakFinding -d ${TplDir}/R0090_FOD_Warp.nii.gz -v -p R0090 -N 0 -n 2 -o ${TplOutputDir}
fi      

if [[ ! -f ${OutputDir}/${ID}0.nii.gz ]]; then 
      peakFinding -d ${InputDir}/${ID}_FOD_Warp.nii.gz -v -p ${ID} -N 0 -n 2 -o ${OutputDir}
fi

if [[ ! -f ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz && ! -f ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz ]]; then 
      RelabelPeaks.py ${ID} ${OutputDir} ${TplOutputDir}/R00900.nii.gz
else 
   echo "One or both relabled peak files exist"      
fi

#Compute GC Scalars( first and second peak ).
if [[ ! -f ${OutputDir}/${ID}_relabled_1st_peak_gc.nii.gz && -f ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz ]]; then 
      GeodesicConcentration -N 1 -a ${InputDir}/${ID}_FOD_Warp.nii.gz -b ${OutputDir}/${ID}_relabeled_1st_peak.nii.gz -p ${OutputDir}/${ID}_relabeled_1st_peak
fi

if [[ ! -f ${OutputDir}/${ID}_relabed_2nd_peak_gc.nii.gz && -f ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz ]]; then 
      GeodesicConcentration -N 1 -a ${InputDir}/${ID}_FOD_Warp.nii.gz -b ${OutputDir}/${ID}_relabeled_2nd_peak.nii.gz -p ${OutputDir}/${ID}_relabeled_2nd_peak
fi

#Compute GFA
if [[ ! -f ${OutputDir}/${ID}_gfa.nii.gz ]]; then 
      echo "Computing GFA for subject " ${ID}
      GeneralFractionalAnisotropy -a ${InputDir}/${ID}_FOD_Warp.nii.gz -p ${OutputDir}/${ID} -N 1
fi
