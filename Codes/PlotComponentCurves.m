function yhat = PlotComponentCurves(x, y, t, c, parameter)
try
	fontSize = 20;
	% Get the means and widths.
	means = parameter(1 : 2 : end);
	widths = parameter(2 : 2 : end);
	% Now plot results.
	hFig2 = figure;
	hFig2.Name = 'Fitted Component Curves';
% 	plot(x, y, '--', 'LineWidth', 2)
	hold on;
	yhat = zeros(1, length(t));
	numGaussians = length(c);
	legendStrings = cell(numGaussians + 2, 1);
	for k = 1 : numGaussians
		% Get each component curve.
		thisEstimatedCurve = c(k) .* gaussian(t, means(k), widths(k));
		% Plot component curves.
		plot(x, thisEstimatedCurve, '-', 'LineWidth', 2);
		hold on;
		% Overall curve estimate is the sum of the component curves.
		yhat = yhat + thisEstimatedCurve;
		legendStrings{k} = sprintf('Estimated Gaussian %d', k);
	end
	% Plot original summation curve, that is the actual curve.
	plot(x, y, 'r-', 'LineWidth', 1)
	% Plot estimated summation curve, that is the estimate of the curve.
	plot(x, yhat, 'k--', 'LineWidth', 2)
	grid on;
	xlabel('X', 'FontSize', fontSize)
	ylabel('Y', 'FontSize', fontSize)
	caption = sprintf('Estimation of %d Gaussian Curves that will fit data.', numGaussians);
	title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
	grid on
	legendStrings{numGaussians+1} = sprintf('Actual original signal');
	legendStrings{numGaussians+2} = sprintf('Sum of all %d Gaussians', numGaussians);
	legend(legendStrings);
	xlim(sort([x(1) x(end)]));
	hFig2.WindowState = 'maximized';
	drawnow;
	
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
end % of PlotComponentCurves