%% Simulation 1:

% Parameters:

Ks = 100; %lbs/in
Kt = 1000; % lbs/in
Ws = 1000; % lbs
Wu = 100; % lbs
Ct = 0.01; % lbs-sec/in

g = 386.06; % in/sec^2

ms = Ws / g; % Sprung mass lbs*sec^2/in
mu = Wu / g; % Unsprung mass in lbs*sec^2/in

freq_low_bound = 4*2*pi; % 4 Hz = 25.13 rad/sec
freq_up_bound = 8*2*pi; % 8 Hz = 50.27 rad / sec

Cs_array = logspace(log10(0.1), log10(100), 10); % In lbs/in/sec

% Legend arrays:
legendEntries = cell(1, 2*length(Cs_array) + 2); % Initialize the cell array
legendCounter = 1; % Initialize a counter for legend entries

legendEntriesBottom = cell(1, length(Cs_array) + 2); 
legendCounterBottom = 1; % Initialize a counter for legend entries

figure(1); % Create the first figure outside the loop
subplot(2, 1, 1)
xscale log;
grid on;
hold on;

subplot(2, 1, 2)
xscale log;
grid on;
hold on;

for i=1:length(Cs_array)

    Cs = Cs_array(i);

    % plot:
    subplot(2, 1, 1);
    [mgs, mgu, mg_iso, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu);
    semilogx(w, mgs, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;

    semilogx(w, mgu, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendEntries{legendCounter} = sprintf('Unsprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;

    % plot:
    subplot(2, 1, 2);
    mg_iso_db = 20*log(mg_iso);
    semilogx(w, mg_iso_db, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendEntriesBottom{legendCounterBottom} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounterBottom = legendCounterBottom + 1;

end

% Initial guess for Cs
Cs_initial = 3.0; % Starting point for the search

options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); % FIXME: Use me

A = []; b = []; Aeq = []; Beq = []; lb = 1.0; ub = 100;
[Cs_optimal, cost_optimal] = fmincon(@objectiveFunction, Cs_initial, A, b, Aeq, Beq, lb, ub, [], options); % TODO: Add options to show optimization by iter

% Output the optimal Cs
fprintf('Optimal Cs: %f, with a cost of %f\n', Cs_optimal, cost_optimal);

subplot(2, 1, 1);
[mgs, mgu, mg_iso, w] = twomass_rel_damp(Ks, Kt, Cs_optimal, Ct, ms, mu);
mg_iso_db = 20*log(mg_iso);

semilogx(w, mgs, 'o-', LineWidth=3.0); % Plot the magnitude vs frequency for the sprung mass
legendEntries{legendCounter} = sprintf('Sprung mass transm. for selected opt. damper Cs = %.2f, Ct = %.2f', Cs_optimal, Ct);
legendCounter = legendCounter + 1;

semilogx(w, mgu, '--', LineWidth=3.0); % Plot the magnitude vs frequency for the unsprung mass
legendEntries{legendCounter} = sprintf('Unsprung mass transm. for selected opt. damper Cs = %.2f, Ct = %.2f', Cs_optimal, Ct);
legendCounter = legendCounter + 1;
title('Relative damping transmissibility plot varying Cs');
legend(legendEntries{1:2*length(Cs_array) + 2}); % Create the legend for the first figure
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;

subplot(2, 1, 2);
semilogx(w, mg_iso_db, '--', LineWidth=3.0); % Plot the magnitude vs frequency for the unsprung mass
legendEntriesBottom{legendCounterBottom} = sprintf('Sprung mass isolation func. for selected opt. damper Cs = %.2f, Ct = %.2f', Cs_optimal, Ct);
legendCounterBottom = legendCounterBottom + 1;

% human_sensitivity = 4*ones(length(w), 1);
% semilogx(w, human_sensitivity, '--', LineWidth=3);
legendEntriesBottom{legendCounterBottom} = sprintf('Human sensitivity');
legend(legendEntriesBottom{:}); % Create the legend for the first figure
xline(freq_low_bound);
xline(freq_up_bound);
ylabel('magnitude (dB)')
hold off;


%% The system with the optimal damper is loaded as a state-space in Matlab for
% Bode plots:

% Mathematical model:

Cs = Cs_optimal;

A_qc = [ 
     0,  1,  0,  0;
    -Ks/ms, -Cs/ms,  Ks/ms,  Cs/ms;
     0,   0,   0,   1;
     Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  -(Cs+Ct)/mu];
 
B_qc = [
    0, 0;
    0, 0;
    0, 0;
    Kt/mu, Ct/mu];
 
C_qc = [
    1, 0, 0, 0;
    0, 0, 1, 0];
 
D_qc = [
    0, 0;0, 0];

% Assuming the optimal damper value and other parameters are defined
% Generate the state-space model
sys_qc = ss(A_qc, B_qc(:,1), C_qc(2,:), D_qc(1,1));

figure;
hold on;
bode(sys_qc);
grid on;
hold off;

%% Isolation function:

C_isolation = [0, 1, 0, 0];

% Bode:
figure;
hold on;
xscale log;
grid on;
[mag, phase, wout] = bode(sys_qc);
magdB = 20*log10(mag); % Convert magnitude to dB

semilogx(squeeze(wout), squeeze(magdB), 'o-', LineWidth=3.0); % Plot the magnitude vs frequency for the sprung mass
yline(4);
xline(freq_low_bound);
xline(freq_up_bound);


%% Set up time domain simulation:
t_initial = 0; t_final = 10;
dt = 0.001; % 1000 Hz
steps = (t_final - t_initial) / dt;
time_space = linspace(t_initial, t_final, steps);
x_initial = zeros(4, 1);

% Parameters:
params = struct('Ks', Ks, 'Kt', Kt, 'Ws', Ws, 'Wu', Wu, 'Ct', Ct, 'Cs', Cs_optimal, 'g', g, 'ms', ms, 'mu', mu);


%% Create an input that has a 1 inch step at 1 second:
amplitude = 1; % inch

road_input1 = amplitude*ones(steps, 1);

road_input1(1:1/dt) = 0;

% Start loop:
x_array = qc_time_simulation(t_initial, x_initial, steps, dt, road_input1, params);

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
x_array = qc_time_simulation(t_initial, x_initial, steps, dt, road_input2, params);


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
x_array = qc_time_simulation(t_initial, x_initial, steps, dt, road_input3, params);


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

amplitude = 0.25; % inch
freq = 10; % Hz
to_rad_sec = 2*pi;

road_input4 = amplitude*sin(freq*to_rad_sec*time_space);

road_input4(1:1/dt) = 0;

% Start loop:
x_array = qc_time_simulation(t_initial, x_initial, steps, dt, road_input4, params);


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
x_array = qc_time_simulation(t_initial, x_initial, steps, dt, road_input5, params);


% Plot x_s:
plot(time_space, x_array(1, :), LineWidth=3.0); 
hold on;
plot(time_space, road_input5, LineWidth=2.0);
legend(["Sprung mass displacement (in)", "Road input (in)"]);
grid on;
xlabel("Time (seconds)")
ylabel("Displacement (in)")
hold off;



%% Functions:

function [x_array] = qc_time_simulation(t0, x0, steps, dt, road_input, params)
    
    Ks = params.Ks;
    ms = params.ms;
    Cs = params.Cs;
    mu = params.mu;
    Kt = params.Kt;
    Ct = params.Ct;

    t = t0; 
    x_array = zeros(4, steps); % 4 rows, steps columns to store the time simulation states over time
    x = x0;

    A = [ 
         0,  1,  0,  0;
        -Ks/ms, -Cs/ms,  Ks/ms,  Cs/ms;
         0,   0,   0,   1;
         Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  -(Cs+Ct)/mu];
     
    B = [
        0, 0;
        0, 0;
        0, 0;
        Kt/mu, Ct/mu];


    for iter = 1:steps
        dx = A*x + B(:,1)*road_input(iter);
        x_array(:, iter) = x;
        x = x + dx*dt;
        t = t + dt;
    end

end


% Objective function that the optimization algorithm will minimize
function cost = objectiveFunction(Cs)
    % Parameters:
    
    Ks = 100; %lbs/in
    Kt = 1000; % lbs/in
    Ws = 1000; % lbs
    Wu = 100; % lbs
    Ct = 0.01; % lbs-sec/in
    
    g = 386.06; % in/sec^2
    
    ms = Ws / g; % Sprung mass lbs*sec^2/in
    mu = Wu / g; % Unsprung mass in lbs*sec^2/in

    % Assume twomass_rel_damp is modified to accept Cs and return a cost metric
    [mgs, mgu, mg_iso, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu);
    mg_iso_db = 20*log(mg_iso);
    % Cost function formulation:
    relevantIndices = (w >= 4*2*pi) & (w <= 8*2*pi); 
    relevantIndices2 = (w >= 4) & (w <= 8);
    cost = 2*mean(mgs(relevantIndices)) + 10*mean(mgs(relevantIndices2) + 0.075*mean(abs(mg_iso(relevantIndices))));
    % cost = mean(mg_iso(relevantIndices)); %max(mgs(relevantIndices));Cs = Cs_array(i);
end


function [mgs, mgu, mg_iso, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu)

    % Mathematical model:
    A=[ 0,  1,  0,  0;
        -Ks/ms, -Cs/ms,  Ks/ms,  Cs/ms;
         0,   0,   0,   1;
         Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  -(Cs+Ct)/mu];
     
    B=[0, 0;
        0, 0;
        0, 0;
        Kt/mu, Ct/mu];
     
    C=[1, 0, 0, 0;
        0, 0, 1, 0];
     
    D=[0, 0;0, 0];

    %twomass calculates the frequency response of a two-mass
    [mag, phase, w]=bode(ss(A,B(:,1), C(1,:), [0]),logspace(0,3)); % Outputs Xs
     mgs(1:50)=mag;
     
    [mag, phase, w]=bode(ss(A,B(:,1), C(2,:), [0]),logspace(0,3)); % Outputs Xu
    mgu(1:50)=mag;

    [mag, phase, w]=bode(ss(A,B(:,1), [0, 1, 0, 0], [0]),logspace(0,3)); % Outputs Xu
    mg_iso(1:50)=mag;
end

