#! /bin/sh

if [ "X$1" != "X" ] ; then
ID=$1
fi

if [ "X$ID" = "X"  ] ; then
echo "No ID given" 1>&2
exit 1
fi

echo " Subject: " ${ID}

#InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}/reg2R0090/normFODReg
InputDir=/sbia/home/raoarvin/comp_space/GeoVariance_Maps/${ID}/Norm_SA_ODF_Reg 
OutputDir=/sbia/home/raoarvin/comp_space/GeoVariance_Maps/${ID}/Norm_SA_ODF_Reg 

#if [[ ! -d ${OutputDir} ]]; then 
#    echo "creating directory for subject " ${ID}
#    mkdir ${OutputDir}
#    mkdir ${OutputDir}/Norm_FOD
#fi

if [[ ! -f ${OutputDir}/${ID}_Norm_SA_mult_peak_char.hdr ]]; then 

echo "Calculating Peak Map Char of unsmoothed FOD"
PeakMapsMakerTest ${InputDir}/${ID}_sa_odf.hdr ${OutputDir}/${ID}_unsmoothed_SA
DiscreteGfaTest ${InputDir}/${ID}_sa_odf.hdr ${OutputDir}/${ID}_unsmoothed_SA
fi

# if [[ -f ${OutputDir}/${ID}_Norm_FOD_peak_maps.hdr && ! -f ${OutputDir}/${ID}_Norm_FOD_first_peak_GC.hdr ]]; then 

#    echo "Spliting the just calculated Peak Maps File"
#    SplitVectorImage.py ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD_peak_maps.hdr ${OutputDir}/${ID}_40it_c2_smoothed_Norm_FOD
#  fi
