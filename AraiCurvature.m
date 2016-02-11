function [parameters]=AraiCurvature(x,y)
%
% Function for calculating the radius of the best fit circle to a set of 
% x-y coordinates.
% Paterson, G. A., (2011), A simple test for the presence of multidomain
% behaviour during paleointensity experiments, J. Geophys. Res., in press,
% doi: 10.1029/2011JB008369
%
% parameters(1) = k
% parameters(2) = a
% parameters(3) = b
% parameters(4) = SSE (goodness of fit)

% Reshape vectors for suitable input
x=reshape(x, length(x), 1);
y=reshape(y, length(y), 1);

% Normalizevectors
x=x./max(x);
y=y./max(y);

% Provide the intitial estimate
E1=TaubinSVD([x,y]);

% Determine the iterative solution
E2=LMA([x,y], E1);

estimates=[E2(3), E2(1), E2(2)];

% Define the function to be minimized and calculate the SSE
func=@(v) sum((sqrt((x-v(2)).^2+(y-v(3)).^2)-v(1)).^2);
SSE=func(estimates);

if (E2(1)<=mean(x) && E2(2)<=mean(y))
    k=-1/E2(3);
else
    k=1/E2(3);
end

parameters=[k; E2(1); E2(2); SSE];

return