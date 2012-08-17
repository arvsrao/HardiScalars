

%
%   [stopTime minR pp W] = CutStreamLine( line )
%
%   line = original streamline points from mrtrix
%   stopTime = time point the curve passes the cutting plan 
%   pp = cut curve, parametrized [0,1]
%   W = sample points of curve
%
function [pp W T] = CutStreamLine( line )

	results=[];
	N=100000;
	timeStep=.00001;
	c=[25.12 31.85 2.06]; %center of seed ball
	d=[23.61 27.98 1.94]; %point on surface ball after which most streamline points have started.

       
        %bspline approximation
	T = linspace(0,1,length(line));
	breaks = linspace(0,1,25);
	pp=splinefit(T,line', breaks);
	
	for ii=1:N
	    results(ii) = abs((d-c)*(ppval(pp,ii*timeStep)'-c)');
	end

        minR=min(results);
        stopTime = find(results<=minR)/N;
        clear results;
        
        %sample the line starting at the point it passes the cutting plane
%        sampleTimes=linspace(stopTime,1,300);
%        W=ppval(pp,sampleTimes);
%        T=linspace(0,1,300);
%        pp = splinefit(T,W,breaks);
        
        %cut the ends off
        c = [11.73 4.76 1.94]; % center of include ball
        d = [13.53 8.18 1.94]; % point on surface of include ball, by which point most streamline have passed.
        results=[];
        
        for ii=1:N
	     results(ii) = abs((d-c)*(ppval(pp, ii*timeStep)'-d)');
        end
	
	minR=min(results);
        endTime = find(results<=minR)/N;
        clear results;
        
	%sample the line starting at the point it passes the cutting plane
        sampleTimes=linspace(stopTime,endTime,500);
        W=ppval(pp,sampleTimes);
        T=linspace(0,1,500);
        pp = splinefit(T,W,breaks);
        
        
        
