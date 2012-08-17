#! /bin/sh

if [ "X$1" != "X" ] ; then
       ID=$1
fi

if [ "X$ID" = "X"  ] ; then
    echo "No ID given" 1>&2
    exit 1
fi

InputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}
OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}

if [[ -f ${InputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr && "$ID" != "R0090" ]]; then
 	
    # if [[ ! -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr  && -f ${OutputDir}/${ID}_FOD_RshPowers.hdr ]]; then 
         
          echo "Smoothing FOD "
          #Accident: many, if not all, of the 20it images were overwritten with 30it infomation, but we're recomputing them 
	#  FODSmooth ${InputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr ${OutputDir}/${ID}_FOD_RshPowers.hdr ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr 10 1
    # fi

     if [[ -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr && ! -f ${OutputDir}/${ID}_smoothed_FOD_30it_first_peak_mag.hdr ]]; then 
	     
         echo "Calculating Peak Map of smoothed FOD"
         PeakImageMaker ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr ${OutputDir}/${ID}_smoothed_FOD_30it
     fi
    
     if [[ -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr && ! -f ${OutputDir}/${ID}_smoothed_FOD_30it_first_peak_GC.hdr ]]; then 
	     
         echo "Calculating GC Map of smoothed FOD"
         GeoVarianceFirstPeakWithOutMasking ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr ${OutputDir}/${ID}_smoothed_FOD_30it
     fi
     
     if [[ -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr && ! -f ${OutputDir}/${ID}_smoothed_FOD_30it_GFA.hdr ]]; then     
         echo "Calculating GFA Map of smoothed FOD"
         DiscreteGfaTest ${OutputDir}/${ID}_FOD_2mm_jLMMSE_warped_RAS_20it_c1.hdr ${OutputDir}/${ID}_smoothed_FOD_30it
     fi       
fi 

if [[ -f ${InputDir}/${ID}_FOD_2mm_jLMMSE_RAS_20it_c1.hdr && "$ID" == "R0090" ]]; then
 	
     if [[ ! -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr  && -f ${OutputDir}/${ID}_FOD_RshPowers.hdr ]]; then 
         
          echo "Smoothing template FOD "
	  FODSmooth ${InputDir}/${ID}_FOD_2mm_jLMMSE_RAS_20it_c1.hdr ${OutputDir}/${ID}_FOD_RshPowers.hdr ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr 10 1
     fi

     if [[ -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr && ! -f ${OutputDir}/${ID}_smoothed_FOD_30it_first_peak_mag.hdr ]]; then 
	     
         echo "Calculating Peak Map of smoothed FOD template"
         PeakImageMaker ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr ${OutputDir}/${ID}_smoothed_FOD_30it
     fi
    
     if [[ -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr && ! -f ${OutputDir}/${ID}_smoothed_FOD_30it_first_peak_GC.hdr ]]; then 
	     
         echo "Calculating GC Map of smoothed FOD template"
         GeoVarianceFirstPeakWithOutMasking ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr ${OutputDir}/${ID}_smoothed_FOD_30it
     fi
     
     if [[ -f ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr && ! -f ${OutputDir}/${ID}_smoothed_FOD_30it_GFA.hdr ]]; then     
         echo "Calculating GFA Map of smoothed FOD template"
          DiscreteGfaTest ${OutputDir}/${ID}_FOD_2mm_jLMMSE_RAS_30it_c1.hdr ${OutputDir}/${ID}_smoothed_FOD_30it
     fi       
fi
