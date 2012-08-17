#! bash/bin

# Make 2nd Peak mask from Template FOD 2nd Peak Map.
# 
#
#


# Compute the peak magnitudes from the 2nd peak image of R0090, the template. 
3dcalc -a 'R0090_relabeled_2nd_peak.nii.gz[0]' -b 'R0090_relabeled_2nd_peak.nii.gz[1' -c 'R0090_relabeled_2nd_peak.nii.gz[2]' -expr 'sqrt(a*a +b*b+c*c)' -prefix '2nd_pm.nii.gz'

# compute intra-image mean and variance for template 2nd peak magnitude.
# 'var' and 'mean' are the sample mean and variance respectively. 
3dBrickStat -non-zero -mean 2nd_pm.nii.gz 
3dBrickStat -non-zero -var 2nd_pm.nii.gz 

#compute the t-stat image. 
3dcalc -a '2nd_pm.nii.gz' -expr '(a-0.137694)/0.00636424' -prefix '2nd_pm_tstat.nii.gz'

#threshold 2nd peak image based on t-stat, at the 68% level. 
3dcalc -a  '2nd_pm_tstat.nii.gz' -expr 'step(a-1)' -prefix '2nd_peak_mask.nii.gz'
