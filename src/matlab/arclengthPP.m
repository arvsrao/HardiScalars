%
%  Usage:  y = arclengthPP(T, pp, dt)
%  where:
%
%  T = 1-d vector of sample points
%  pp = polyline representing the line.
%  dt = time step.

function y = arclengthPP(T, pp, dt)

   
   W=ppval(pp,T);
   qq=ppdiff(pp,1);
   dW=ppval(qq,T);
   y = dt*sum(sqrt(diag(dW'*dW)));
   

%  How to make interpolated line.  
%
%  breaks = linspace(0,1,20);
%  T = linspace(0,1,length(lineOne));
%  pp=splinefit(T,lineOne, breaks);


