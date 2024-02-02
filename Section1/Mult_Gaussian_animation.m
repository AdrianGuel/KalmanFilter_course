% Parameters
mu = [0, 0];         % Mean of the Gaussian distribution
sigma = [1, 0.5;     % Covariance matrix
         0.5, 1];
num_points = 1e4;    % Number of points to sample
animation_speed = 0.01; % Speed of animation (adjust as needed)

% Generate animation frames
fig = figure;
for i = 1:num_points
    % Sample from the multivariate Gaussian distribution
    U = chol(sigma);
    point = mu+U'*randn(2,1);
    %point = mvnrnd(mu, sigma);
    
    % Plot the point
    plot(point(1), point(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    title(['Frame ' num2str(i) '/' num2str(num_points)]);
    xlabel('X-axis');
    ylabel('Y-axis');
    axis equal;
    grid on;
    hold on;
    
    % Pause for animation
    pause(animation_speed);
    
    % Clear previous point to create animation effect
    %if i < num_points
    %    clf(fig);
    %end
end
