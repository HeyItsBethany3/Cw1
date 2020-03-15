% Reads in data
value = csvread('q.csv');

% Retrieves q

num = value(1,1);
qValue = value(2,1:end-1);
fValue = value(3,1:end-1);


% Analtical f
f = @(x) (exp((x.^2)+x).*sin(pi.*(5/4).*x));

% Create nodes

xNodes = zeros(1,num+1);
h = 2/num;
for i=0:num
    xNodes(i+1) = (-1)+(i*h);
end

% Create plot
figure(1);
plot(xNodes, fValue, 'c', xNodes, qValue, 'm');
xlabel('x'); legend('f(x)', 'p(x)');
title("Least squares approximation p(x) to f(x)");
ax = gca; ax.FontSize = 14; axis tight;
