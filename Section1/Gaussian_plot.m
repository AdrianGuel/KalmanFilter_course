% Parameters for the Gaussian distribution
mu = 0;         % Mean
sigma = 1;      % Standard deviation

% Generate data points for the Gaussian distribution
x = linspace(-5, 5, 1000); % Adjust the range as needed
y = (1/(sigma * sqrt(2 * pi))) * exp(-(x - mu).^2 / (2 * sigma^2));

% Generate random points using randn
numPoints = 1000; % Adjust the number of points as needed
randomPoints = mu + sigma * randn(1, numPoints);

% Plot the Gaussian distribution
plot(x, y, 'LineWidth', 2);

hold on; % Add the following line to hold the current plot

% Plot the random points using randn
plot(randomPoints, normpdf(randomPoints, mu, sigma), 'r.', 'MarkerSize', 10);

hold off; % Release the hold on the plot

% Add labels and title
xlabel('X-axis');
ylabel('Probability Density');
title('Gaussian Distribution with Random Points');

% Add a grid to the plot
grid on;

% Display a legend
legend('Gaussian Distribution', 'Random Points');
