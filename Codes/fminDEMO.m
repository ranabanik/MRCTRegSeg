optionsDemo = optimset('PlotFcns',@optimplotfval); 
% optimset-> matlab function
% 'PlotFcns' -> function variable
% @optimplotfval -> just an input to the varible to see the plots
baikka = @(xDemo)100*(xDemo(2) - xDemo(1)^2)^2 + (1 - xDemo(1))^2;
% Rosenbrock's equation: 100*(x2 - x1^2)^2 + (1 - x1)^2;
% here the value of variables x1 and x2 will change to find the best fit. 
% the values are written in a single variable xDemo = [x1, x2]
% the function handle @(xDemo) means to change the xDemo inside
% Rosenbrock's equation. 
x0 = [-1.2,1];
% this is the starting point. 
xDemo = fminsearch(baikka,x0,options);