% Parameters
numSteps = 100;   % Number of steps
deltaT = 0.1;     % Time step size
sigma = 1;        % Standard deviation

% Initialize position arrays
xPosition = zeros(1, numSteps);
yPosition = zeros(1, numSteps);

% Create a figure
figure;
xlabel('X Position');
ylabel('Y Position');
title('2D Brownian Motion');

% Initialize the plot with the starting point
hPlot = plot(0, 0, 'LineWidth', 2);
hold on;
grid on;

for step = 2:numSteps
    % Generate random increments using normal distribution for x and y
    incrementX = sigma * sqrt(deltaT) * randn;
    incrementY = sigma * sqrt(deltaT) * randn;
    
    % Update the positions
    xPosition(step) = xPosition(step-1) + incrementX;
    yPosition(step) = yPosition(step-1) + incrementY;
    
    % Update the plot
    set(hPlot, 'XData', xPosition(1:step), 'YData', yPosition(1:step));
    
    % Pause to control the animation speed
    pause(0.1);
    
    % Force the plot to update
    drawnow;
end

hold off;
