% Reads in data
value = csvread('q.csv');

% Retrieves q
qValue = value(1,1:end-1);

f = @(x) (exp((x.^2)+x)*sin(pi*1.25*x));

% Create nodes
num = 101;
xNodes = zeros(num);
h = 2.0/num;
for i=0:100
    xNodes(i+1) = (-1)+(i*h);
end

% Creates plot
figure(1);
plot(xNodes, f(xNodes), 'c', xNodes, qValue, 'm');
