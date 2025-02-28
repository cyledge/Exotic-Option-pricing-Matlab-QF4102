% Test script for the FD_ids_amCall function

% Given parameters
S0 = 1.231;            % Current price of the underlying asset
X = 1;                 % Strike price
r = 0.03;              % Risk-free interest rate (3%)
sigma = 0.2;           % Volatility (20%)
q = 0.04;              % Continuous dividend yield (4%)
T = 3/12;              % Time to maturity in years (3 months)
xmin = -4;             % Minimum value of the grid
xmax = 4;              % Maximum value of the grid
N = 1000;              % Number of time steps
omega = 1.3;           % Relaxation parameter for SOR
eps = 1e-6;            % Error tolerance

% Test cases for different numbers of spatial grid points (I)
I_values = 80:80:800;  % Array of grid point counts
%N = 3;
%I_values = [80, 160];

% Preallocate result array
option_values = zeros(size(I_values));

% Run the function for each value of I and store the results
for idx = 1:length(I_values)
    I = I_values(idx);
    option_values(idx) = FD_ids_amCall(S0, X, r, T, sigma, q, I, N, xmin, xmax, omega, eps);
    fprintf('Option value for I = %d: %f\n', I, option_values(idx));
end

% Plot the estimated option values versus I
figure;
plot(I_values, option_values, '-o');
xlabel('Number of spatial grid points (I)');
ylabel('Estimated American Call Option Value');
title('Estimated Option Value vs Number of Grid Points');
grid on;