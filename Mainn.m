% Main script for 2D Ising Model Monte Carlo Simulation
clear; clc; close all;
tic

% Simulation parameters
L = 100;                % Lattice size
J = 1;                  % Coupling constant
kB = 1;                 % Boltzmann constant
Tc = 2/log(1+sqrt(2));  % Critical temperature (~2.269)
n_thermal = 2e4;        % Thermalization sweeps
n_eq = 1e4;             % Equilibration sweeps
n_avg = 1e4;            % Averaging sweeps

% Initialize all figures
for figNum = 1:8
    figure(figNum);
    clf;
    set(gcf, 'Position', [100+50*figNum, 100+50*figNum, 600, 500]);
end

% ========================
% Thermalization Analysis (Figures 1-2)
% ========================
fprintf('Running thermalization analysis...\n');
T_thermal = [2.0, 2.5];  % Below and above Tc
IC_options = {'all_up', 'all_down', 'random'};
colors = ['r', 'g', 'b'];

for T_idx = 1:length(T_thermal)
    T = T_thermal(T_idx);
    
    for ic = 1:length(IC_options)
        % Initialize lattice
        spin = initializeLattice(L, IC_options{ic});
        
        % Run Metropolis with measurement
        [spin, Energy, Magnetization] = Metropolis(spin, T, L, n_thermal, true);
        
        % Plot thermalization
        figure(1);
        subplot(2,1,T_idx);
        plot(Magnetization, colors(ic), 'DisplayName', IC_options{ic});
        title(sprintf('Magnetization (T=%.1f)', T));
        xlabel('Sweeps'); ylabel('magnetisation'); grid on;
        hold on;
        
        figure(2);
        subplot(2,1,T_idx);
        plot(Energy, colors(ic), 'DisplayName', IC_options{ic});
        title(sprintf('Energy (T=%.1f)', T));
        xlabel('Sweeps'); ylabel('e(t)'); grid on;
        hold on;
    end
end

% Add legends
for fig = 2:1
    figure(fig);
    for subp = 2:1
        subplot(2,1,subp);
        legend('show', 'Location', 'best');
    end
end

% ========================
% Microscopic Configurations (Section 1.3)
% ========================
fprintf('Generating microscopic configurations...\n');
T_config = [2.0, 2.5];  % Below and above Tc
snap_times = 2e3:2e3:18e3; % Exactly 9 snapshots at 2000, 4000,...,18000 sweeps

for T_idx = 1:length(T_config)
    T = T_config(T_idx);
    spin = initializeLattice(L, 'random'); % Start from random configuration
    
    % Create figure for this temperature
    fig_num = 3 + (T == 2.5); % Figure 3 for T=2.0, Figure 4 for T=2.5
    figure(fig_num); clf;
    sgtitle(sprintf('Microscopic Configurations at T=%.1f', T));
    
    % Run simulation and capture snapshots
    snapshot_count = 1;
    for sweep = 1:18e3 % Run for 18000 sweeps
        spin = Metropolis(spin, T, L, 1, false); % Do 1 sweep at a time
        
        % Capture snapshot at specified times
        if ismember(sweep, snap_times)
            subplot(3,3,snapshot_count);
            imagesc(spin);
            colormap([1 1 0; 0 0 0]); % White=+1, Black=-1
            title(sprintf('t = %d sweeps', sweep));
            set(gca,'XTick',[],'YTick',[]);
            axis square;
            snapshot_count = snapshot_count + 1;
        end
    end
end

% ========================
% Phase Transition Analysis (Figures 5-8)
% ========================
fprintf('Running phase transition analysis...\n');
T_range = linspace(1, 4, 31);  % Temperature range

% Preallocate arrays
Energy_avg = zeros(size(T_range));
Magnetization_avg = zeros(size(T_range));
MagSus = zeros(size(T_range));
SpecHeat = zeros(size(T_range));

% Start from ordered configuration
spin = ones(L, L);

for i = 1:length(T_range)
    T = T_range(i);
    fprintf('Processing T = %.2f (%d/%d)\n', T, i, length(T_range));
    
    % Equilibration (no measurement)
    spin = Metropolis(spin, T, L, n_eq, false);
    
    % Measurement (with data collection)
    [spin, Energy, Magnetization] = Metropolis(spin, T, L, n_avg, true);
    
    % Calculate averages
    Energy_avg(i) = mean(Energy);
    Magnetization_avg(i) = mean(abs(Magnetization));
    MagSus(i) = (L^2/T)*(mean(Magnetization.^2) - mean(Magnetization)^2);
    SpecHeat(i) = (L^2/T^2)*(mean(Energy.^2) - mean(Energy)^2);
end

% Plot phase transition results
plotPhaseTransitionResults(T_range, Energy_avg, Magnetization_avg, MagSus, SpecHeat, Tc);

toc

% ========================
% Helper Functions
% ========================
function spin = initializeLattice(L, IC_type)
    switch IC_type
        case 'all_up'
            spin = ones(L, L);
        case 'all_down'
            spin = -ones(L, L);
        case 'random'
            spin = sign(rand(L, L) - 0.5);
    end
end

function plotPhaseTransitionResults(T_range, Energy_avg, Magnetization_avg, MagSus, SpecHeat, Tc)
    % Figure 5: Magnetization
    figure(5);
    plot(T_range, Magnetization_avg, 'b-o', 'LineWidth', 1.5);
    hold on; grid on;
    plot([Tc Tc], [0 1], 'r--', 'LineWidth', 1.5);
    title('Mean Magnetization vs Temperature');
    xlabel('Temperature'); ylabel('Magnetisation');
    legend('MC Data', 'T_c', 'Location', 'best');
    
    % Figure 6: Energy
    figure(6);
    plot(T_range, Energy_avg, 'b-o', 'LineWidth', 1.5);
    hold on; grid on;
    plot([Tc Tc], ylim, 'r--', 'LineWidth', 1.5);
    title('Mean Energy vs Temperature');
    xlabel('Temperature'); ylabel('Energy');
    legend('MC Data', 'T_c', 'Location', 'best');
    
    % Figure 7: Specific Heat
    figure(7);
    plot(T_range, SpecHeat, 'm-*', 'LineWidth', 1.5);
    hold on; grid on;
    plot([Tc Tc], ylim, 'r--', 'LineWidth', 1.5);
    title('Specific Heat vs Temperature');
    xlabel('Temperature'); ylabel('Spesific Heat,c');
    legend('MC Data', 'T_c', 'Location', 'best');
    
    % Figure 8: Susceptibility
    figure(8);
    plot(T_range, MagSus, 'm-*', 'LineWidth', 1.5);
    hold on; grid on;
    plot([Tc Tc], ylim, 'r--', 'LineWidth', 1.5);
    title('Magnetic Susceptibility vs Temperature');
    xlabel('Temperature'); ylabel('Magnetic Susceptibility');
    legend('MC Data', 'T_c', 'Location', 'best');
end