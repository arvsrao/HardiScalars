#! /bin/sh

if [ "X$1" != "X" ] ; then
       ID=$1
fi

if [ "X$ID" = "X"  ] ; then
    echo "No ID given" 1>&2
    exit 1
fi

#normalized FOD images

InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}/reg2R0090/normFODReg
#InputDir=/sbia/sbiaprj/autism/TimRoberts/Processed_HARDI_RAS/${ID}/reg2R0090/normFODReg #Tim's dataset
#OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}

#if [[ ! -d ${OutputDir} ]]; then 
#    echo "creating directory for subject " ${ID}
#    mkdir ${OutputDir}
#    mkdir ${OutputDir}/Norm_FOD
#fi

OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}/Norm_FOD


if [[ -f ${InputDir}/${ID}_Norm_FOD_Warp.nii.gz ]]; then
 	 
    # if [[ ! -f ${OutputDir}/${ID}_FOD_RshPowers.nii.gz ]]; then  
        
	#  echo "Calculating FOD Powers Image"
        #  computeRshPowers -d ${InputDir}/${ID}_Norm_FOD_Warp.nii.gz -p ${ID}_FOD_RshPowers -o ${OutputDir}	  
    # fi

     if [[ -f ${OutputDir}/${ID}_FOD_RshPowers.nii.gz && ! -f ${OutputDir}/${ID}_Norm_FOD_Warp_40it_c2.hdr ]]; then 
        
          echo "Smoothing FOD"
	  FODSmooth ${InputDir}/${ID}_Norm_FOD_Warp.nii.gz ${OutputDir}/${ID}_FOD_RshPowers.nii.gz ${OutputDir}/${ID}_Norm_FOD_Warp_40it_c2.hdr 40 2
     fi

     if [[ -f ${OutputDir}/${ID}_Norm_FOD_Warp_40it_c2.hdr && ! -f ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD_peak_maps.hdr ]]; then 
	     
         echo "Calculating Peak Map of smoothed FOD"
         PeakMapsMaker ${OutputDir}/${ID}_Norm_FOD_Warp_40it_c2.hdr ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD
         DiscreteGfaTest ${OutputDir}/${ID}_Norm_FOD_Warp_40it_c2.hdr ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD 
     fi
     
     if [[ -f ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD_peak_maps.hdr && ! -f ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD_first_peak_GC.hdr ]]; then 
	     
         echo "Spliting the just calculated Peak Maps File"
         SplitVectorImage.py ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD_peak_maps.hdr ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD
     fi
fi
