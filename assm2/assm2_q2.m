%% for 2(iii)

S0 = 100;         % Initial stock price
X = 100;          % Strike price
r = 0.03;         % Risk-free rate
T = 1;            % Time to maturity
sigma = 0.22;     % Volatility of the asset
q = 0;            % No dividends
N = 4;            % Number of time steps
L = 2;            % Discretization parameter

% Call the function with these parameters
option_price = fsg_fixArithAsianCallNew(S0, X, r, T, sigma, q, N, L);

% Display the result
disp(['The price of the fixed-strike arithmetic Asian call option is: ', num2str(option_price)]);


%% for 2(v)

% Define parameters for the arithmetic average Asian call option
S0 = 95;           % Current underlying price
X = 90;            % Strike price
r = 0.04;          % Risk-free rate
T = 0.5;           % Time to expiry (in years)
sigma = 0.30;      % Volatility of the underlying
q = 0;             % Dividend yield
N_values = [60, 120, 180, 240];  % Number of time periods
rho_values = [1, 0.5, 0.25];     % Values for rho
runavg = 93;       % Historical average of the underlier

% Initialize storage for results
results = zeros(length(N_values), length(rho_values));
runtimes = zeros(length(N_values), length(rho_values));

% Loop over values of N and rho to compute option value estimates
for n_idx = 1:length(N_values)
    N = N_values(n_idx);
    
    for rho_idx = 1:length(rho_values)
        L = 1 / rho_values(rho_idx);
        
        % Start timer to record runtime
        tic;
        
        % Call the function to get the option value
        Nhist = 0;  % Set Nhist for current calculation
        v = fsg_fixArithAsianCall(S0, X, r, T, sigma, q, N, L, runavg, Nhist);
        
        % Stop timer and record runtime
        runtimes(n_idx, rho_idx) = toc;
        
        % Store the result
        results(n_idx, rho_idx) = v;
    end
end

% Display the results
fprintf('Option Value Estimates:\n');
disp(array2table(results, 'VariableNames', {'Rho_1', 'Rho_0_5', 'Rho_0_25'}, 'RowNames', {'N_60', 'N_120', 'N_180', 'N_240'}));

fprintf('\nRuntimes (seconds):\n');
disp(array2table(runtimes, 'VariableNames', {'Rho_1', 'Rho_0_5', 'Rho_0_25'}, 'RowNames', {'N_60', 'N_120', 'N_180', 'N_240'}));


%% for 2(vi)

% Script to plot the runtimes versus N for the function fsg_fixArithAsianCall

% Define parameters for the arithmetic average Asian call option
S0 = 95;           % Current underlying price
X = 90;            % Strike price
r = 0.04;          % Risk-free rate
T = 0.5;           % Time to expiry (in years)
sigma = 0.30;      % Volatility of the underlying
q = 0;             % Dividend yield
N_values = [60, 120, 180, 240];  % Number of time periods
rho_values = [1, 0.5, 0.25];     % Values for rho
runavg = 93;       % Historical average of the underlier
Nhist = 2;         % Number of historical time periods used to calculate runavg

% Initialize storage for runtimes
runtimes = zeros(length(N_values), length(rho_values));

% Loop over values of N and rho to compute runtimes
for n_idx = 1:length(N_values)
    N = N_values(n_idx);
    
    for rho_idx = 1:length(rho_values)
        L = round(1 / rho_values(rho_idx));  % Ensure L is an integer
        
        % Start timer to record runtime
        tic;
        
        % Call the function to get the option value (not storing results)
        try
            v = fsg_fixArithAsianCall(S0, X, r, T, sigma, q, N, L, runavg, Nhist);
        end
        
        % Stop timer and record runtime
        runtimes(n_idx, rho_idx) = toc;
    end
end

% Plotting the runtimes versus N
figure;
hold on;
for rho_idx = 1:length(rho_values)
    plot(N_values, runtimes(:, rho_idx), '-o', 'DisplayName', sprintf('\rho = %.2f', rho_values(rho_idx)));
end
hold off;

xlabel('Number of Time Periods (N)');
ylabel('Runtime (seconds)');
title('Runtime vs Number of Time Periods for Different \rho Values');
legend;
grid on;

