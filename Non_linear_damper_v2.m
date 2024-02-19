%% Sim1 Non-linear damper (different rebound and compression coefficients):

% Parameters:

Ks = 100; %lbs/in
Kt = 1000; % lbs/in
Ws = 1000; % lbs
Wu = 100; % lbs
Ct = 0.01; % lbs-sec/in

g = 386.06; % in/sec^2

ms = Ws / g; % Sprung mass lbs*sec^2/in
mu = Wu / g; % Unsprung mass in lbs*sec^2/in


%% Set up time domain simulation:
t_initial = 0; t_final = 10;
dt = 0.0005; % 100 Hz
steps = (t_final - t_initial) / dt;
time_space = linspace(t_initial, t_final, steps);
Cs_optimal = 50;
x_initial = zeros(4, 1);

% Parameters:
params = struct('Ks', Ks, 'Kt', Kt, 'Ws', Ws, 'Wu', Wu, 'Ct', Ct, 'Cs', Cs_optimal, 'g', g, 'ms', ms, 'mu', mu);

%% Create an input that has a 1 inch step at 1 second:
amplitude = 1; % inch

road_input1 = amplitude*ones(steps, 1);

road_input1(1:1/dt) = 0;

% Start loop:
x_array = qc_time_nonlinear_simulation(t_initial, x_initial, steps, dt, road_input1, params);

% Plot x_s:
figure;
plot(time_space, x_array(1, :), LineWidth=3.0); 
hold on;
plot(time_space, road_input1, LineWidth=2.0);
legend(["Sprung mass displacement (in)", "Road input (in)"]);
grid on;
xlabel("Time (seconds)")
ylabel("Displacement (in)")
hold off;

%% Create an input that has 1/4 inch step at 1 second:
amplitude = 1/4; % inch

road_input2 = amplitude*ones(steps, 1);

road_input2(1:1/dt) = 0;

% Start loop:
x_array = qc_time_nonlinear_simulation(t_initial, x_initial, steps, dt, road_input2, params);


% Plot x_s:
figure;
plot(time_space, x_array(1, :), LineWidth=3.0); 
hold on;
plot(time_space, road_input2, LineWidth=2.0);
legend(["Sprung mass displacement (in)", "Road input (in)"]);
grid on;
xlabel("Time (seconds)")
ylabel("Displacement (in)")
hold off;

%% Create an input +- sine wave of amplitude 1 at 1 Hz:

amplitude = 1; % inch
freq = 1; % Hz
to_rad_sec = 2*pi;

road_input3 = amplitude*sin(freq*to_rad_sec*time_space);

road_input3(1:1/dt) = 0;

% Start loop:
x_array = qc_time_nonlinear_simulation(t_initial, x_initial, steps, dt, road_input3, params);


% Plot x_s:
figure;
plot(time_space, x_array(1, :), LineWidth=3.0); 
hold on;
plot(time_space, road_input3, LineWidth=2.0);
legend(["Sprung mass displacement (in)", "Road input (in)"]);
grid on;
xlabel("Time (seconds)")
ylabel("Displacement (in)")
hold off;

%% Create an input +- sine wave of amplitude 1/4 at 10 Hz:

amplitude = 1; % inch
freq = 10; % Hz
to_rad_sec = 2*pi;

road_input4 = amplitude*sin(freq*to_rad_sec*time_space);

road_input4(1:1/dt) = 0;

% Start loop:
x_array = qc_time_nonlinear_simulation(t_initial, x_initial, steps, dt, road_input4, params);


% Plot x_s:
figure;
plot(time_space, x_array(1, :), LineWidth=3.0); 
hold on;
plot(time_space, road_input4, LineWidth=2.0);
legend(["Sprung mass displacement (in)", "Road input (in)"]);
grid on;
xlabel("Time (seconds)")
ylabel("Displacement (in)")
hold off;

%% A 2 inch pothole 3 ft wide at 50 mph

width = 3*12; % in
amplitude = 2; % in
speed = 50; % mph
speed_in_sec = speed*17.6; % in/sec
step_time_width = width / speed_in_sec;
step_time = 1;

road_input5 = zeros(steps, 1);
road_input5(step_time/dt:(step_time + step_time_width)/dt) = -2; % -2 inches of displacement due to the 2 in pothole

% Start loop:
x_array = qc_time_nonlinear_simulation(t_initial, x_initial, steps, dt, road_input5, params);


% Plot x_s:
figure;
plot(time_space, x_array(1, :), LineWidth=3.0); 
hold on;
plot(time_space, road_input5, LineWidth=2.0);
legend(["Sprung mass displacement (in)", "Road input (in)"]);
grid on;
xlabel("Time (seconds)")
ylabel("Displacement (in)")
hold off;

function [x_array] = qc_time_nonlinear_simulation(t0, x0, steps, dt, road_input, params)
    
    Ks = params.Ks;
    ms = params.ms;
    Cs_compression = 1.25*params.Cs; % Assuming 25% more damping in compression
    Cs_rebound = params.Cs; % Standard damping in rebound
    mu = params.mu;
    Kt = params.Kt;
    Ct = params.Ct;

    t = t0; 
    x_array = zeros(4, steps); % 4 rows, steps columns to store the time simulation states over time
    x = x0;
    for iter = 1:steps

        spring_deflection = x(3) - x(1);  
        if spring_deflection > 3.0
           spring_deflection = 3.0; 
        end
        if spring_deflection < -3.0
            spring_deflection = -3.0;
        end

        % Modify damping force based on compression or rebound
        velocity_difference = x(4) - x(2);
        if velocity_difference > 0
            % Compression
            Fc = Cs_compression * velocity_difference;
        else
            % Rebound
            Fc = Cs_rebound * velocity_difference; % Using negative velocity difference directly
        end

        ddx_s = (1/ms)*(Ks*(spring_deflection) + Fc);
        ddx_u = (1/mu)*(Kt*(road_input(iter) - x(3)) - Ks*(spring_deflection) - Fc);
        x(2) = x(2) + dt*ddx_s;
        x(4) = x(4) + dt*ddx_u;
        x(1) = x(1) + dt*x(2);
        x(3) = x(3) + dt*x(4);
        x_array(:, iter) = x;
        t = t + dt;

    end

end