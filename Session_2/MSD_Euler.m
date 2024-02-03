% Mass-Spring-Damper System Parameters
m = 1;          % Mass (kg)
k = 10;         % Spring constant (N/m)
c = 0.5;        % Damping coefficient (Ns/m)

% Simulation Parameters
dt = 0.01;      % Time step (s)
t_final = 5;    % Final simulation time (s)

% Initial Conditions
x0 = 1;         % Initial position (m)
v0 = 0;         % Initial velocity (m/s)

% Initialize arrays for time, position, and velocity
t = 0:dt:t_final;
x = zeros(size(t));
v = zeros(size(t));

% Initial conditions
x(1) = x0;
v(1) = v0;

% Euler method simulation
for i = 2:length(t)
    % Compute acceleration (from Newton's second law)
    a = (-c*v(i-1) - k*x(i-1))/m;
    
    % Update velocity and position using Euler method
    v(i) = v(i-1) + a*dt;
    x(i) = x(i-1) + v(i-1)*dt;
end

% Plot the results
figure;

subplot(2,1,1);
plot(t, x, 'b-', 'LineWidth', 1.5);
title('Mass-Spring-Damper System');
xlabel('Time (s)');
ylabel('Position (m)');
grid on;

subplot(2,1,2);
plot(t, v, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

sgtitle('Position and Velocity vs Time');

