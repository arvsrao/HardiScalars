#! /bin/sh

if [ "X$1" != "X" ] ; then
	ID=$1
fi

if [ "X$ID" = "X"  ] ; then
	echo "No ID given" 1>&2
	exit 1
fi

InputDir=/sbia/comp_space/raoarvin/GeoVariance_Maps/${ID}
FixedImage=/sbia/comp_space/raoarvin/GeoVariance_Maps/R0090

if [ -f ${InputDir}/${ID}_wm_mask.hdr ]; then 
      echo ${ID} " mask flirted with flirt matrix"
else 
      echo "flirt mask with R0090 mask"
      flirt -in ${InputDir}/${ID}_T1_Reg_seg_2.hdr -ref ${FixedImage}/R0090_T1_Reg_seg_2.hdr -init ${InputDir}/${ID}_flirt_mat -out ${InputDir}/${ID}_wm_mask -applyxfm -interp nearestneighbour 
fi

if [ "$ID" != "R0090" ]; then
      echo "applying Deformation field from ODF registration."
      applyDeformationField -d ${InputDir}/${ID}_wm_mask.hdr -f ${InputDir}/${ID}_odf_reg_defField.hdr -o ${InputDir} -p ${ID}_wm_mask_warped -i 1
fi


