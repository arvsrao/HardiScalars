%
%  Compute Max Cluster Size for Random Relabeling.
%
%      USAGE:  maxClusterSize = MaxClusterSizeStat(Y, Ymask, vol, group, age) 
%  
%
function maxClusterSize = MaxClusterSizeStat(Y, Ymask, vol, group, age)

	N=length(group);
	
	%Permute Labels
	labelPerm=randperm(N);
	AgePerm = term(age(labelPerm));
	GroupPerm = term( var2fac( group(labelPerm), {'patient'; 'control'}) );
	contrast=GroupPerm.control-GroupPerm.patient;
	
	%Compute t-stats
	slm = SurfStatLinMod(Y, 1 + GroupPerm + AgePerm);
	slm = SurfStatT(slm, contrast);

	%Compute p-values
	tau=57/116;
	pMap=1-fcdf(slm.t.*slm.t*tau, 2, 57);

	tmapMask=(pMap<=.05);
	map=reshape(tmapMask.*Ymask, [128,128,80]);

	CC = bwconncomp(map);
	stats = regionprops(CC, 'Area');

	areas=[];
	for ii=1:length(stats)
	    areas = [stats(ii).Area areas];
	end

	maxClusterSize = max(areas);  % max cluster size statistic 
	
