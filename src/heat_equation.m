% Heat Equation Simulation 
% Author: Justin-Marian
% Date: 2024
% Description: This script simulates the 1D heat equation for metal bars.
% Considered metals: Aluminum, Iron, and Copper.
% The simulation discretizes space and time for stability.
% Neumann boundary conditions are assumed.
%
% MIT License
%
% Copyright 2024 Justin-Marian
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
% IN THE SOFTWARE.
%
function heat_equation()
    % Simulates the 1D heat equation for metal bars.
    % Considered metals: Aluminum, Iron, and Copper.
    %
    % Output:
    %  - The simulation discretizes space and time for stability.
    %  - Neumann boundary conditions are assumed.
    %  - Designed for educational purposes.
    %
    clc;
    close all;

    %% Physical parameters of the metal bars

    % Length in meters of the metal bars (Al, Fe, Cu)
    L1 = 0.04;
    L2 = 0.07;
    L3 = 0.05;

    % Thermal conductivity in W/mK
    lambda1 = 237; % Al
    lambda2 = 80;  % Fe
    lambda3 = 401; % Cu

    % Density in kg/m^3
    ro1 = 2700; % Al
    ro2 = 7860; % Fe
    ro3 = 8960; % Cu

    % Specific heat in J/kgK
    c1 = 900; % Al
    c2 = 250; % Fe
    c3 = 380; % Cu

    % Calculate thermal diffusivity for each metal bar
    D1 = lambda1 / ro1 / c1;
    D2 = lambda2 / ro2 / c2;
    D3 = lambda3 / ro3 / c3;

    % Half lengths of the metal bars
    L1_half = L1 / 2;
    L2_half = L2 / 2;
    L3_half = L3 / 2;
    
    % Total length of the metal bars
    L_total = L1_half + L2_half + L3_half;
    
    % Discretization of the metal bars
    P = 111;
    
    %%

    % Start from the left end of the second bar
    x = linspace(-L1 - L2_half, L_total - L3_half, P);
    dx = x(2) - x(1);

    % Characteristic time
    Tc1 = L1^2 / D1;
    Tc2 = L2^2 / D2;
    Tc3 = L3^2 / D3;
    Tc = mean([Tc1, Tc2, Tc3]);

    % Discretization of time
    t0 = 0; tf = 5 * Tc;
    N = 10^6; t = linspace(t0, tf, N);
    dt = t(2) - t(1);

    % Initial temperatures and diffusivity
    T1 = 220; % Kelvin (Al)
    T2 = 77;  % Temperature (Fe)
    T3 = 340; % (Cu)

    %%

    % Vector of initial temperatures
    T0 = zeros(1, P);
    index_L1_end = round(L1 / (L1 + L2 + L3) * P);
    index_L2_end = round((L1 + L2) / (L1 + L2 + L3) * P);

    T0(1:index_L1_end) = T1;
    T0(index_L1_end + 1:index_L2_end) = T2;
    T0(index_L2_end + 1:end) = T3;
    
    % Set diffusion coefficient based on the material type
    D = D1 * ones(1, P);
    D(index_L1_end + 1:index_L2_end) = D2;
    D(index_L2_end + 1:end) = D3;

    % Matrix to store temperatures at each point and time
    T = zeros(N, P);
    T(1, :) = T0;

    % Calculate the temperature distribution over time
    for i = 1:N - 1
        for j = 2:P - 1
            T(i + 1, j) = T(i, j) + dt / dx^2 *...
                (D(j + 1) * (T(i, j + 1) -...
                 T(i, j)) -...
                 D(j) * (T(i, j) - T(i, j - 1)));
        end
        
        % Neumann boundary conditions
        T(i + 1, 1) = T(i + 1, 2);
        T(i + 1, P) = T(i + 1, P - 1);
    end

    %%

    % Dynamic simulation
    figure('Name', 'Heat Equation Simulation');
    simt = 0;

    while simt <= tf
        % Find the nearest time index to the current time
        index = min(abs(t - simt)) == abs(t - simt);
    
        % Clear previous plot and update the temperature distribution
        hold off; %
        plot(x(1:index_L1_end), T(index, 1:index_L1_end), 'r--o', 'DisplayName', 'Al');
        hold on;

        plot(x(index_L1_end+1:index_L2_end), T(index, index_L1_end+1:index_L2_end), 'g--o', 'DisplayName', 'Fe');
        plot(x(index_L2_end+1:end), T(index, index_L2_end+1:end), 'b--o', 'DisplayName', 'Cu');

        % Add lines indicating the boundaries between materials
        line([-L1, -L1], [T2, T3], 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'L1 <=> L2');
        line([L1_half / 2, L1_half / 2], [T2, T3], 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'L2 <=> L3');

        xlabel('x/m');
        ylabel('T/K');
        axis([-L1 - L2_half, L_total - L3_half, T2, T3]);
    
        % Update position and string value of the text
        text(x(45), 1.25 * T(index, 45), ['t = ', num2str(round(t(index) * 10)), ' ds'], 'Color', [0.9290 0.6940 0.1250], 'FontSize', 12);

        % Add legend
        legend('Location', 'southeast');
    
        % Pause for a short time to allow visualization
        pause(1e-3);
        simt = simt + 0.1;
    end
end
