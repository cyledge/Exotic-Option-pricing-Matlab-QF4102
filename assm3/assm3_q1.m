% Script to test both MC_assm3 and MC_assm3_CV functions

% Parameters for the test
S0 = [12, 10];        % Initial prices of the two assets
X = 12;               % Strike price
r = 0.04;             % Risk-free rate (4%)
T = 2;                % Time to maturity (2 years)
sigma = [0.25, 0.35]; % Volatilities of the two assets
q = [0.02, 0.03];     % Dividend yields of the two assets
rho = -0.12;          % Correlation coefficient
H = 14;               % Barrier level
N = 500;              % Number of time periods

% Simulation settings for different sample sizes
P_values = [300, 3000, 30000]; % Number of sample paths for estimation
num_runs = 30;                 % Number of independent simulation runs

% Preallocate result arrays for both methods
results_MC = zeros(length(P_values), num_runs);    % For MC_assm3
results_MC_CV = zeros(length(P_values), num_runs); % For MC_assm3_CV

% Loop through each P value and perform 30 independent runs for MC_assm3
for p_idx = 1:length(P_values)
    P = P_values(p_idx);
    for run = 1:num_runs
        % Call the MC_assm3 function
        results_MC(p_idx, run) = MC_assm3(S0, X, r, T, sigma, rho, q, N, H, P);
    end
end

% Loop through each P value and perform 30 independent runs for MC_assm3_CV
for p_idx = 1:length(P_values)
    P = P_values(p_idx);
    for run = 1:num_runs
        % Call the MC_assm3_CV function
        results_MC_CV(p_idx, run) = MC_assm3_CV(S0, X, r, T, sigma, rho, q, N, H, P);
    end
end

% Calculate the mean and standard error for each P value for both methods
mean_estimates_MC = mean(results_MC, 2);
standard_errors_MC = std(results_MC, 0, 2) / sqrt(num_runs);

mean_estimates_MC_CV = mean(results_MC_CV, 2);
standard_errors_MC_CV = std(results_MC_CV, 0, 2) / sqrt(num_runs);

% Display the results in a table for both methods
table_MC = table(P_values', mean_estimates_MC, standard_errors_MC, ...
    'VariableNames', {'P', 'Mean_Estimate_MC', 'Standard_Error_MC'});

table_MC_CV = table(P_values', mean_estimates_MC_CV, standard_errors_MC_CV, ...
    'VariableNames', {'P', 'Mean_Estimate_MC_CV', 'Standard_Error_MC_CV'});

% Display the tables
disp('Results of Monte Carlo Estimation (Standard Method):');
disp(table_MC);

disp('Results of Monte Carlo Estimation with Control Variate:');
disp(table_MC_CV);
