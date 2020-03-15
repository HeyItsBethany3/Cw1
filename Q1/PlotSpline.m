% Reads in data
value = csvread('PlotSpline.csv');

% Retrieves h and error values
hvec = value(1,1:end-1);
error = value(2, 1:end-1);

% Plots error against h values
figure(1);
loglog(hvec, error, 'm');
xlabel('h'); ylabel('Error');
title('Spline error as step-size h changes for point x=1/3');

% Finds gradient of plot
coefficients = polyfit(log(hvec), log(error), 1);
slope = coefficients(1);
disp(slope); % Approximately 4
