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

freq_low_bound = 4*2*pi; % 4 Hz = 25.13 rad / sec
freq_up_bound = 8*2*pi; % 8 Hz = 50.27 rad / sec

Cs_array = logspace(log10(0.1), log10(100), 10); % In lbs/in/sec

% Legend arrays:
legendEntriesTop = cell(1, 2*length(Cs_array) + 2); % Initialize the cell array
legendCounterTop = 1; % Initialize a counter for legend entries

legendEntriesBottom = cell(1, length(Cs_array) + 1); % Initialize the cell array
legendCounterBottom = 1; % Initialize a counter for legend entries


figure(1); % Create the first figure outside the loop

subplot(2,1,1);
xscale log;
grid on;
hold on;

subplot(2,1,2);
xscale log;
grid on;
hold on;

for i=1:length(Cs_array)

    Cs = Cs_array(i);

    % plot:
    subplot(2,1,1);
    [mag_sprung, mag_unsprung, mag_isolation, w1, w2] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu);
    semilogx(w1, mag_sprung, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendEntriesTop{legendCounterTop} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounterTop = legendCounterTop + 1;

    semilogx(w1, mag_unsprung, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendEntriesTop{legendCounterTop} = sprintf('Unsprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounterTop = legendCounterTop + 1;

    % plot:
    subplot(2,1,2);
    semilogx(w2, mag_isolation, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendEntriesBottom{legendCounterBottom} = sprintf('Sprung mass accel Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounterBottom = legendCounterBottom + 1;

end

% Initial guess for Cs
Cs_initial = 40.0; % Starting point for the search

options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); % FIXME: Use me

A = []; b = []; Aeq = []; Beq = []; lb = 1.0; ub = 100;
[Cs_optimal, cost_optimal] = fmincon(@objectiveFunction, Cs_initial, A, b, Aeq, Beq, lb, ub); % TODO: Add options to show optimization by iter

% Print the optimal Cs
fprintf('Optimal Cs: %f, with a cost of %f\n', Cs_optimal, cost_optimal);

% Generate frequency responses:
[mag_sprung, mag_unsprung, mag_isolation, w1, w2] = twomass_rel_damp(Ks, Kt, Cs_optimal, Ct, ms, mu);

subplot(2, 1, 1); % Top plot
semilogx(w1, mag_sprung, 'o-', LineWidth=3.0); % Plot the magnitude vs frequency for the sprung mass
legendEntriesTop{legendCounterTop} = sprintf('Sprung mass transm. for selected opt. damper Cs = %.2f, Ct = %.2f', Cs_optimal, Ct);
legendCounterTop = legendCounterTop + 1;

semilogx(w1, mag_unsprung, '--', LineWidth=3.0); % Plot the magnitude vs frequency for the unsprung mass
legendEntriesTop{legendCounterTop} = sprintf('Unsprung mass transm. for selected opt. damper Cs = %.2f, Ct = %.2f', Cs_optimal, Ct);
legendCounterTop = legendCounterTop + 1;

title('Relative damping transmissibility plot varying Cs');
legend(legendEntriesTop{1:2*length(Cs_array) + 2}); % Create the legend for the first figure
% semilogx(w, human_sensitivity, '--', LineWidth=2)
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;


subplot(2, 1, 2); % Bottom plot
semilogx(w2, mag_isolation, 'o-', LineWidth=3.0); % Plot the magnitude vs frequency for the sprung mass
legendEntriesBottom{legendCounterBottom} = sprintf('Sprung mass isolation func. for selected opt. damper Cs = %.2f, Ct = %.2f', Cs_optimal, Ct);
legendCounterBottom = legendCounterBottom + 1;

title('Relative damping isolation function plot varying Cs');
legend(legendEntriesBottom{1:length(Cs_array) + 1}); % Create the legend for the first figure
% semilogx(w, human_sensitivity, '--', LineWidth=2)
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;


% human_sensitivity = power(10,4/20)*ones(length(w), 1);


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
figure;
sys_qc = ss(A_qc, B_qc(:,1), C_qc(1,:), D_qc(1,1));
bode(sys_qc);
% figure;
% [mag, phase, wout] = bode(sys_qc);
% magdB = 20*log10(mag); % Convert magnitude to dB
% 
% h = bodeplot(sys_qc);

%% Isolation function:

C_isolation = [0, 1, 0, 0];

% Bode:
figure;
hold on;
xscale log;
[mag, phase, wout] = bode(sys_qc);
magdB = 20*log10(mag); % Convert magnitude to dB

semilogx(squeeze(wout), squeeze(magdB), 'o-', LineWidth=3.0); % Plot the magnitude vs frequency for the sprung mass
yline(4);
xline(freq_low_bound);
xline(freq_up_bound);


%% Time domain simulation:
t_initial = 0; t_final = 10;
dt = 0.01; % 100 Hz
steps = (t_final - t_initial) / dt;
time_space = linspace(t_initial, t_final, steps);

% Create an input that has one inch step at 0.5 Hz:
amplitude = 1; % inch
freq_hz = 0.5;
road_freq = freq_hz*2*pi;
period = 1 / road_freq;
road_input = amplitude*sin(road_freq*time_space).*ones(1, steps);

% zero_indices = road_input < 0;
% one_indices = road_input >= 0;

% road_input(zero_indices) = 0;
% road_input(one_indices) = 1;

x_initial = zeros(4,1);
x = zeros(4,1);

% Parameters:
params = struct('Ks', Ks, 'Kt', Kt, 'Ws', Ws, 'Wu', Wu, 'Ct', Ct, 'Cs', Cs_optimal, 'g', g, 'ms', ms, 'mu', mu);

% Start loop:
x_array = qc_time_simulation(t_initial, x_initial, steps, dt, road_input, params);

% Plot x_s:
plot(time_space, x_array(1, :)); 
hold on;
plot(time_space, road_input);
hold off;

%% Non-linear damper selection:




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
    [mgs_sus, mgu_sus, mag_iso, w1, w2] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu);
    
    % Cost function formulation:
    relevantIndices = (w1 >= 4*2*pi) & (w1 < 8*2*pi); % rad/sec for the 4 to 8 Hz region % TODO: Might need to add human sensitivity here
    relevantIndices2 = (w1 >= 4) & (w1 < 8); % (rad/sec) More or less within this range is where the ride natural frequency occurs
    cost = mean(mgs_sus(relevantIndices)) + mean(mgu_sus(relevantIndices2)); %max(mgs(relevantIndices));
end


function [mgs_sus, mgu_sus, mgs_iso, w1, w2] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu)

    % Mathematical model:
    A=[ 0,  1,  0,  0;
        -Ks/ms, -Cs/ms,  Ks/ms,  Cs/ms;
         0,   0,   0,   1;
         Ks/mu,  Cs/mu,  -(Ks+Kt)/mu,  -(Cs+Ct)/mu];
     
    B=[0, 0;
        0, 0;
        0, 0;
        Kt/mu, Ct/mu];
     
    % C=[1, 0, 0, 0;
    %     0, 0, 1, 0];

    C = eye(4);
     
    D=[0, 0; 0, 0];

    % Generate frequency responses:

    % Bode from road displacement to sprung mass displacement:
    [mag_sprung_sus, phase, w1] = bode(ss(A, B(:,1), [1, 0, 0, 0], D(1,1)), logspace(0,3)); % Outputs Xs
    mgs_sus(1:50) = mag_sprung_sus;

    % Bode from road displacement to usprung mass displacement:
    [mag_unsprung_sus, phase, w1] = bode(ss(A, B(:,1), [0, 0, 1, 0], D(1,1)), logspace(0,3)); % Outputs Xu
    mgu_sus(1:50) = mag_unsprung_sus;
    
    % Bode from road input velocity to sprung mass acceleration (isolation function):
    [mag_iso, phase, w2] = bode(ss(A, B(:,1), [0, 1, 0, 0], D(1,1)), logspace(0,3)); % Outputs Xu
    mgs_iso(1:50) = mag_iso;
end

