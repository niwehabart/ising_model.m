function [spin, Energy, Magnetization] = Metropolis(spin, T, L, n_sweeps, do_measure)
% Metropolis algorithm for 2D Ising model
% Inputs:
%   spin - initial spin configuration
%   T - temperature
%   L - lattice size
%   n_sweeps - number of Monte Carlo sweeps
%   do_measure - whether to record measurements
% Outputs:
%   spin - final spin configuration
%   Energy - energy time series (if do_measure)
%   Magnetization - magnetization time series (if do_measure)

% Initialize outputs
if do_measure
    Energy = zeros(1, n_sweeps);
    Magnetization = zeros(1, n_sweeps);
else
    Energy = [];
    Magnetization = [];
end

for sweep = 1:n_sweeps
    % One sweep (L^2 spin flip attempts)
    for k = 1:L^2
        % Random spin selection
        row = randi(L);
        col = randi(L);
        
        % Periodic boundary conditions
        above = mod(row - 2, L) + 1;
        below = mod(row, L) + 1;
        left = mod(col - 2, L) + 1;
        right = mod(col, L) + 1;
        
        % Calculate energy change
        dE = 2 * spin(row,col) * (spin(above,col) + spin(row,right) + ...
                                  spin(row,left) + spin(below,col));
        
        % Metropolis acceptance
        if dE <= 0 || rand() < exp(-dE/T)
            spin(row,col) = -spin(row,col);
        end
    end
    
    % Measure observables if requested
    if do_measure
        % Calculate energy
        energy = 0;
        for x = 1:L
            for y = 1:L
                above = mod(x - 2, L) + 1;
                right = mod(y, L) + 1;
                energy = energy - spin(x,y)*(spin(above,y) + spin(x,right));
            end
        end
        
        % Calculate magnetization
        magnetization = sum(spin, 'all');
        
        % Store measurements
        Energy(sweep) = energy / L^2;
        Magnetization(sweep) = magnetization / L^2;
    end
end
end