function MaxClusterSizeStat_core(count, infile, OUTPATH)

	load(infile);
        N=length(group);
	reset(RandStream.getDefaultStream,sum(100*clock));
	
	%Permute Labels
	labelPerm=randperm(N);
	AgePerm = term(age(labelPerm));
	GenderPerm = term(gender(labelPerm));
	GroupPerm = term( var2fac( group(labelPerm), {'patient'; 'control'}) );
	contrast=GroupPerm.control-GroupPerm.patient;
	
	%Compute t-statsx
	slm = SurfStatLinMod(Y, 1 + GroupPerm + AgePerm + GenderPerm);
	slm = SurfStatT(slm, contrast);
        clear Y 
        
	%Compute p-values
	tau=57/116;
	pMap=1-fcdf(slm.t.*slm.t*tau, 2, 57);

	tmapMask=(pMap<=.05);
	map=reshape(tmapMask.*Ymask, [128,128,80]);
        clear tmapMask Ymask  
         
	CC = bwconncomp(map);
	stats = regionprops(CC, 'Area');
        clear map CC 
        
	areas=[];
	for ii=1:length(stats)
	    areas = [stats(ii).Area areas];
	end

	maxClusterSize = max(areas);  % max cluster size statistic 
	     
	outFile = [ OUTPATH 'TMaxStat_' count '.mat' ];	
	save(outFile, 'maxClusterSize');
       
end
