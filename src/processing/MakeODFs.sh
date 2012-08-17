#! /bin/sh

#Reconstruction Odf from DWI

if [ "X$1" != "X" ] ; then
	ID=$1
fi

if [ "X$ID" = "X"  ] ; then
	echo "No ID given" 1>&2
	exit 1
fi

InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}
InputDir_MPR=/sbia/sbiaprj/autism/TobaccoCAR/Processed/${ID}
OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}

if [ -d $OutputDir ]; then
   echo $OutputDir " already exists."
else 
   echo "Making output directory in sbia/comp_space/raoarvin/GeoVariance_Maps"
   mkdir ${OutputDir}
fi

if [ -f ${OutputDir}/${ID}_B0.hdr ]; then
	echo " B0 already extracted."
else
        echo "extracting B0"
        extractB0.py -d ${InputDir}/${ID}_HARDI-DWI_2mm_2_jLMMSE_merged_RAS.nii -b ${InputDir}/${ID}_HARDI-DWI_2mm_2_jLMMSE_merged_RAS.bval -p ${ID}_B0 -o ${OutputDir}
fi

if [ -f ${OutputDir}/${ID}_MPR_res_correct_RAS.hdr ]; then
        echo " T1 rotated already. "
else
        echo "rotating T1"
        3dresample -orient lpi -prefix ${OutputDir}/${ID}_MPR_res_correct_RAS.hdr -inset ${InputDir_MPR}/${ID}_MPR_res_correct.nii.gz
fi        
        
if [ -f  ${OutputDir}/${ID}_brain_mask.hdr ]; then
        echo " Brain mask already created "
else  
        echo "create whole-brain mask"
        bet ${OutputDir}/${ID}_MPR_res_correct_RAS.hdr ${OutputDir}/${ID}_brain -m 
fi

if [ -f ${OutputDir}/${ID}_T1_Reg.hdr ]; then 
      echo " T1 flirted with B0 already "
else 
      echo "flirting T1 to B0"
      flirt -ref ${OutputDir}/${ID}_B0 -in ${OutputDir}/${ID}_brain -out ${OutputDir}/${ID}_T1_Reg.hdr
fi

if [ -f ${OutputDir}/${ID}_T1_Reg_seg.hdr ]; then 
      echo " T1 already segmented in B0 space. "
else 
      echo "segmenting T1 in B0 space"
      fast --nopve -g -o ${OutputDir}/${ID}_T1_Reg ${OutputDir}/${ID}_T1_Reg.hdr
fi       

if [ -f ${OutputDir}/${ID}_odf.hdr ]; then 
      echo " ODF reconstruction already done."
else 
      echo "Odf Reconstruction"
      odfRecon -d ${InputDir}/${ID}_HARDI-DWI_2mm_2_jLMMSE_merged_RAS.nii -p ${ID}_odf -r ${InputDir}/${ID}_HARDI-DWI_2mm_2_jLMMSE_merged_RAS.bvec -m ${OutputDir}/${ID}_brain_mask.img -l 8 -n 1 -v -o ${OutputDir}
fi

#Now Register images to R0090
InputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}
FixedImage=/sbia/comp_space/raoarvin/GeoVariance_Maps/R0090

if [ -f ${InputDir}/${ID}_flirt_mat ]; then 
      echo ${ID}_B0 " flirted with R0090_B0 already "
else 
      echo "outputing flirt matrix "
      flirt -ref ${FixedImage}/R0090_B0.hdr -in ${InputDir}/${ID}_B0.hdr -omat ${InputDir}/${ID}_flirt_mat
fi

if [ -f ${InputDir}/${ID}_odf_flirt.hdr ]; then
      echo "flirt matrix already applied."
else 
      echo "applying flirt matrix. "
      applyFlirtMatToRshImg -f ${FixedImage}/R0090_odf.img -d ${InputDir}/${ID}_odf.img -m ${InputDir}/${ID}_flirt_mat -p ${ID}_odf_flirt -o ${InputDir}
fi

if [ -f ${InputDir}/${ID}_odf_flirt.hdr ]; then
      echo "odf flirt image exists."
else 
      echo "odf flirt image doesn't exist."
fi

if [ "$ID" != "R0090" ]; then
     echo "Registering " $ID " to R0090" 
     RSHDemonsRegistration -f ${FixedImage}/R0090_odf.img -m ${InputDir}/${ID}_odf_flirt.img -p ${ID}_odf_reg -o ${InputDir} -v
else
     echo "No need to register template (R0090)"
fi
