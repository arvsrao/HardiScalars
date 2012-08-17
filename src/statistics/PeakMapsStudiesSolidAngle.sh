#! /bin/sh

if [ "X$1" != "X" ] ; then
       ID=$1
fi

if [ "X$ID" = "X"  ] ; then
    echo "No ID given" 1>&2
    exit 1
fi

#normalized SA ODF images

InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}
#InputDir=/sbia/sbiaprj/autism/TimRoberts/Processed_HARDI_RAS/${ID}/reg2R0090/normFODReg #Tim's dataset
OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}

if [[ ! -d ${OutputDir} ]]; then 
    echo "creating directory for subject " ${ID}
    mkdir ${OutputDir}
    mkdir ${OutputDir}/Norm_SA_ODF
fi

OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}/Norm_SA_ODF

if [[ -f ${InputDir}/${ID}_Norm_SA_ODF_Warp.nii.gz ]]; then
 	 
     if [[ ! -f ${OutputDir}/${ID}_SA_ODF_RshPowers.nii.gz ]]; then  
        
	  echo "Calculating FOD Powers Image"
          computeRshPowers -d ${InputDir}/${ID}_Norm_SA_ODF_Warp.nii.gz -p ${ID}_SA_ODF_RshPowers -o ${OutputDir}	  
     fi

     if [[ -f ${OutputDir}/${ID}_SA_ODF_RshPowers.nii.gz && ! -f ${OutputDir}/${ID}_Norm_SA_ODF_Warp_30it_c1.hdr ]]; then 
        
          echo "Smoothing FOD "
	  FODSmooth ${InputDir}/${ID}_Norm_SA_ODF_Warp.nii.gz ${OutputDir}/${ID}_SA_ODF_RshPowers.nii.gz ${OutputDir}/${ID}_Norm_SA_ODF_Warp_30it_c1.hdr 30 1
     fi

     if [[ -f ${OutputDir}/${ID}_SA_OFD_Warp_30it_c1.hdr && ! -f ${OutputDir}/${ID}_smoothed_Norm_SA_ODF_peak_maps.hdr ]]; then 
	     
         echo "Calculating Peak Map of smoothed FOD"
         PeakMapsMaker ${OutputDir}/${ID}_Norm_SA_ODF_Warp_30it_c1.hdr ${OutputDir}/${ID}_smoothed_Norm_SA_ODF
     fi
     
     if [[ -f ${OutputDir}/${ID}_smoothed_Norm_SA_ODF_peak_maps.hdr && ! -f ${OutputDir}/${ID}_smoothed_Norm_SA_ODF_first_peak_GC.hdr ]]; then 
	     
         echo "Spliting the just calculated Peak Maps File"
         SplitVectorImage.py ${OutputDir}/${ID}_smoothed_Norm_SA_ODF_peak_maps.hdr ${OutputDir}/${ID}_smoothed_Norm_SA_ODF
     fi
fi
