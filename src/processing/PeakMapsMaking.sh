#! /bin/sh

if [ "X$1" != "X" ] ; then
       ID=$1
fi

if [ "X$ID" = "X"  ] ; then
    echo "No ID given" 1>&2
    exit 1
fi

echo " Subject: " ${ID}
#normalized FOD images

InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}/reg2R0090/normFODReg
#InputDir=/sbia/sbiaprj/autism/TimRoberts/Processed_HARDI_RAS/${ID}/reg2R0090/normFODReg #Tim's dataset
OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}

if [[ ! -d ${OutputDir} ]]; then 
    echo "creating directory for subject " ${ID}
    mkdir ${OutputDir}
    mkdir ${OutputDir}/Norm_FOD
fi

OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}/Norm_FOD
 	 
if [[ ! -f ${OutputDir}/${ID}_unsmoothed_Norm_FOD_peak_maps_char.nii ]]; then 
	     
      echo "Calculating Peak Map Char of unsmoothed FOD"
      GcWithPeakFindingCalculator ${InputDir}/${ID}_Norm_FOD_Warp.nii.gz ${OutputDir}/${ID}_nf_gc_maps.nii.gz
      #DiscreteGfaTest ${InputDir}/${ID}_Norm_FOD_Warp.nii.gz ${OutputDir}/${ID}_unsmoothed_Norm_FOD 
fi
     
#if [[ -f ${OutputDir}/${ID}_Norm_FOD_mult_peak_char.hdr && ! -f ${OutputDir}/${ID}_Norm_FOD_first_peak_GC.hdr ]]; then 
	     
 #     echo "Spliting the just calculated Peak Maps File"
 #     SplitVectorImage.py ${OutputDir}/${ID}_unsmoothed_Norm_FOD_peak_maps.hdr ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD
#fi
