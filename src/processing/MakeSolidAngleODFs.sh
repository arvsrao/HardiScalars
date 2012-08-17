#! /bin/sh

#Reconstruction SA Odf from DWI

if [ "X$1" != "X" ] ; then
	ID=$1
fi

if [ "X$ID" = "X"  ] ; then
	echo "No ID given" 1>&2
	exit 1
fi

InputDir=/sbia/sbiaprj/autism/TobaccoCAR/Processed_HARDI_RAS/${ID}
#InputDir_MPR=/sbia/sbiaprj/autism/TobaccoCAR/Processed/${ID}
OutputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}/Norm_SA_ODF_Reg

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
        3dcalc -prefix ${OutputDir}/${ID}_B0.hdr -a ${InputDir}/${ID}_DWI_merged_RAS.nii.gz[0] -expr 'a'
fi

if [ -f ${OutputDir}/${ID}_T1_res_correct_RAS.hdr ]; then
        echo " T1 rotated already. "
else
        echo "rotating T1"
        3dresample -orient lpi -prefix ${OutputDir}/${ID}_T1_res_correct_RAS.hdr -inset ${InputDir}/${ID}_T1_corrected_RAS.nii.gz
fi        
        
if [ -f  ${OutputDir}/${ID}_brain_mask.hdr ]; then
        echo " Brain mask already created "
else  
        echo "create whole-brain mask"
        bet ${OutputDir}/${ID}_T1_res_correct_RAS.hdr ${OutputDir}/${ID}_brain -m 
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

if [[ -f ${OutputDir}/${ID}_sa_odf.hdr || -f ${OutputDir}/${ID}_sa_odf.hdr.gz || -f ${OutputDir}/${ID}_sa_odf.img.gz ]]; then 
      echo " Solid Angle ODF reconstruction already done."
else 
      echo "Solid Angle Odf Reconstruction"
      diffusionInfo.py -d ${InputDir}/${ID}_DWI_merged_RAS.nii.gz -o ${OutputDir} -s
      
      SolidAngleOdfRecon -d ${InputDir}/${ID}_DWI_merged_RAS.nii.gz -p ${ID}_sa_odf.hdr -r ${OutputDir}/${ID}_DWI_merged_RAS.nii.bvec -m ${OutputDir}/${ID}_brain_mask.img -l 8 -n 0 -v -o ${OutputDir}
fi

#Register later
