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

if [[ -f ${InputDir}/${ID}_flirt_mat || "$ID" = "R0090" ]]; then 
      echo ${ID}_B0 " flirted with R0090_B0 already "
else 
      echo "outputing flirt matrix "
      flirt -ref ${FixedImage}/R0090_B0.hdr -in ${InputDir}/${ID}_B0.hdr -omat ${InputDir}/${ID}_flirt_mat
fi

if [[ -f ${InputDir}/${ID}_odf_flirt.hdr || "$ID" = "R0090" ]]; then
       echo "flirt matrix already applied."
else 
       echo "applying flirt matrix."
      applyFlirtMatToRshImg -f ${FixedImage}/R0090_sa_odf.hdr -d ${InputDir}/${ID}_sa_odf.hdr -m ${InputDir}/${ID}_flirt_mat -p ${ID}_odf_flirt -o ${InputDir}
fi

if [[ -f ${InputDir}/${ID}_odf_flirt.hdr || "$ID" = "R0090" ]]; then
      echo "odf flirt image exists."
else 
      echo "odf flirt image doesn't exist. Registration failed!!"
fi

if [ "$ID" != "R0090" ]; then
     echo "Registering " $ID " to R0090" 
     RSHDemonsRegistration -f ${FixedImage}/R0090_sa_odf.hdr -m ${InputDir}/${ID}_odf_flirt.hdr -p ${ID}_sa_odf_reg -o ${InputDir} -v
else
     echo "No need to register template (R0090)"
fi

if [[ -f ${InputDir}/${ID}_sa_odf_reg.hdr && "$ID" != "R0090" ]]; then
     echo "deleting odf flirt image, 'cause we don't need it anymore" 
     rm ${InputDir}/${ID}_odf_flirt.* -f
else
     echo "registration failed for " ${ID}
fi

if [[ -f ${InputDir}/${ID}_sa_odf_reg.hdr && "$ID" != "R0090" ]]; then
      echo "Zipping sa_odf image and deformation field, 'cause we won't be doing any processing on it" 
      gzip ${InputDir}/${ID}_sa_odf.img ${InputDir}/${ID}_sa_odf.hdr
      gzip ${InputDir}/${ID}_sa_odf_reg_defField.img ${InputDir}/${ID}_sa_odf_reg_defField.hdr
else
      echo "registration failed for " ${ID}
fi

