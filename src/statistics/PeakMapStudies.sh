#! /bin/sh

if [ "X$1" != "X" ] ; then
       ID=$1
fi

if [ "X$ID" = "X"  ] ; then
    echo "No ID given" 1>&2
    exit 1
fi

InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}
OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}

#/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/RXXXX/RXXXX_FOD_2mm_jLMMSE_warped_RAS.nii

if [[ -f ${InputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS.hdr && "$ID" != "R0090" ]]; then
 	
      echo "Calculating Peak Image for SA Odf's"
      PeakImageMaker ${OutputDir}/${ID}_smoothed_sa_odf_10it.hdr ${OutputDir}/${ID}_sa_odf

      echo "Calculating Peak Map of FOD"
      PeakImageMaker ${InputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS.hdr ${OutputDir}/${ID}_FOD
fi

if [[ -f ${InputDir}/R0090_FOD_2mm_jLMMSE_warped_RAS.hdr && "$ID" == "R0090" ]]; then

     echo "Calculating Peak Image for SA Odf Template"
     PeakImageMaker ${OutputDir}/${ID}_smoothed_sa_odf_10it_cond25.hdr ${OutputDir}/${ID}_sa_odf
	      
     echo "Calculating Peak Map of FOD's Template"
     PeakImageMaker ${InputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS.hdr ${OutputDir}/${ID}_FOD

fi
