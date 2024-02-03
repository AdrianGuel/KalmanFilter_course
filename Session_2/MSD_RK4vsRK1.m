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
x_euler = zeros(size(t));
v_euler = zeros(size(t));
x_rk4 = zeros(size(t));
v_rk4 = zeros(size(t));

% Initial conditions
x_euler(1) = x0;
v_euler(1) = v0;
x_rk4(1) = x0;
v_rk4(1) = v0;

% Euler method simulation
for i = 2:length(t)
    a_euler = (-c*v_euler(i-1) - k*x_euler(i-1))/m;  % Acceleration from Newton's second law
    v_euler(i) = v_euler(i-1) + a_euler*dt;          % Update velocity
    x_euler(i) = x_euler(i-1) + v_euler(i-1)*dt;     % Update position
end

% 4th order Runge-Kutta method simulation
for i = 2:length(t)
    % Runge-Kutta update for velocity
    k1v = (-c*v_rk4(i-1) - k*x_rk4(i-1))/m;
    k2v = (-c*(v_rk4(i-1) + 0.5*dt*k1v) - k*(x_rk4(i-1) + 0.5*dt*v_rk4(i-1)))/m;
    k3v = (-c*(v_rk4(i-1) + 0.5*dt*k2v) - k*(x_rk4(i-1) + 0.5*dt*(v_rk4(i-1) + 0.5*dt*k2v)))/m;
    k4v = (-c*(v_rk4(i-1) + dt*k3v) - k*(x_rk4(i-1) + dt*(v_rk4(i-1) + 0.5*dt*k3v)))/m;
    v_rk4(i) = v_rk4(i-1) + (1/6)*(k1v + 2*k2v + 2*k3v + k4v)*dt;
    
    % Runge-Kutta update for position
    k1x = v_rk4(i-1);
    k2x = (v_rk4(i-1) + 0.5*dt*k1x);
    k3x = (v_rk4(i-1) + 0.5*dt*k2x);
    k4x = (v_rk4(i-1) + dt*k3x);
    x_rk4(i) = x_rk4(i-1) + (1/6)*(k1x + 2*k2x + 2*k3x + k4x)*dt;
end

% Plot the results
figure;

subplot(2,1,1);
plot(t, x_euler, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Euler Method');
hold on;
plot(t, x_rk4, 'r-', 'LineWidth', 1.5, 'DisplayName', '4th Order Runge-Kutta');
title('Mass-Spring-Damper System');
xlabel('Time (s)');
ylabel('Position (m)');
grid on;
legend;

subplot(2,1,2);
plot(t, v_euler, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Euler Method');
hold on;
plot(t, v_rk4, 'r-', 'LineWidth', 1.5, 'DisplayName', '4th Order Runge-Kutta');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;
legend;

sgtitle('Position and Velocity vs Time (Comparison)');
