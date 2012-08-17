indir='/home/arvind/Documents/MATLAB/StreamFiles/'


% Creates Plot Streamlines with Colors
 figure;
for ii=1:100
    if exist( [ indir 'scalarval_' num2str(ii) '.txt' ] )~=0
       
       M=dlmread([ indir 'Streamline_bspline_' num2str(ii) '.txt' ] );
       M(:,4)=[];
       color=dlmread([ indir 'scalarval_' num2str(ii) '.txt' ]);
       scatter3(M(:,1),M(:,2), M(:,3),[],color(:,2));
       hold on;
      
    end
end

colorbar;

%Compute Mean & Variance of Streamlines coming from GC values of Streamlines.
colorMean = zeros(500,2);
colorVar = zeros(500,2);
count=0;

for ii=1:100
    if exist( [ indir 'scalarval_' num2str(ii) '.txt' ] )~=0
       count = count +1;       
       color = dlmread([ indir 'scalarval_' num2str(ii) '.txt' ]);
       colorMean = colorMean + color;
       colorVar = colorVar + color.*color;
    end
end

colorMean = colorMean/count; % unbiased sample mean
colorVar = (colorVar - count*colorMean.*colorMean) / (count-1); %unbiased sample variance

%Make ColorMean graph of GC
figure;
T=linspace(0,1,500); %--the x-variable
plot(T', colorMean(:,2));
title('Mean GC Along Streamlines');
ylabel('Mean GC');


%Make ColorMean graph of GFA
figure;
plot(T', colorMean(:,1));
title('Mean GFA Along Streamlines');
ylabel('Mean GFA');


%Make ColorVar graph of GC
figure;
plot(T', colorVar(:,2));
title('Variance of GC Along Streamlines');
ylabel('Variance of GC');


%Make ColorVar graph of GFA
figure;
plot(T', colorVar(:,1));
title('Variance of GFA Along Streamlines');
ylabel('Variance of GFA');


gfaMean = mean(colorMean(:,1))
gcMean = mean(colorMean(:,2))

gfaVar = var(colorMean(:,1))
gcVar = var(colorMean(:,2))












