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
OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}/Norm_FOD

if [[ -f ${OutputDir}/${ID}_FOD_RshPowers.nii.gz && ! -f ${OutputDir}/${ID}_Norm_FOD_Warp_35it_c1.hdr ]]; then 
        
          echo "Smoothing FOD"
	  FODSmooth ${InputDir}/${ID}_Norm_FOD_Warp.nii.gz ${OutputDir}/${ID}_FOD_RshPowers.nii.gz ${OutputDir}/${ID}_Norm_FOD_Warp_35it_c1.hdr 35 1
fi
