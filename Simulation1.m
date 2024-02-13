%% Simulation 1:

% Parameters:

Ks = 100; %lbs/in
Kt = 1000; % lbs/in
Ws = 1000; % lbs
Wu = 100; % lbs
Ct = 0.01; % lbs-sec/in

g = 386.06; % in/sec^2

ms = Ws / g; % Sprung mass lbs*sec^2/in
mu = Ws / g; % Unsprung mass in lbs*sec^2/in


Cs_array = logspace(log10(0.1), log10(100), 10); % In lbs/in/sec

% Legend arrays:
legendEntries = cell(1, 2*length(Cs_array)); % Initialize the cell array
legendCounter = 1; % Initialize a counter for legend entries


figure; % Create the first figure outside the loop
xscale log;
grid on;
hold on;

for i=1:length(Cs_array)

    Cs = Cs_array(i);

    % plot:
    % subplot(2,1,1);
    [mgs, mgu, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu);

    semilogx(w', mgs, 'o-'); % Plot the magnitude vs frequency for the sprung mass
    legendEntries{legendCounter} = sprintf('Sprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;

    semilogx(w, mgu, '--'); % Plot the magnitude vs frequency for the unsprung mass
    legendEntries{legendCounter} = sprintf('Unsprung mass Cs = %.2f, Ct = %.2f', Cs, Ct);
    legendCounter = legendCounter + 1;

end


title('Relative damping transmissibility plot varying Cs');
legend(legendEntries{1:2*length(Cs_array)}); % Create the legend for the first figure
xlabel('frequency [rad/sec]');
ylabel('amplitude ratio');
hold off;


function [mgs, mgu, w] = twomass_rel_damp(Ks, Kt, Cs, Ct, ms, mu)

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
    [mag, phase, w]=bode(ss(A,B(:,1),C(1,:),[0]),logspace(0,3)); % Outputs Xs
     mgs(1:50)=mag;
     
    [mag, phase, w]=bode(ss(A,B(:,1),C(2,:),[0]),logspace(0,3)); % Outputs Xu
    mgu(1:50)=mag;
end


