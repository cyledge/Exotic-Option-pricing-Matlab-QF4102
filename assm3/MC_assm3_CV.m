function v = MC_assm3_CV(S0, X, r, T, sigma, rho, q, N, H, P)
    % Inputs:
    % S0: Vector of initial prices of two assets
    % X: Strike price
    % r: Risk-free interest rate
    % T: Time to maturity (in years)
    % sigma: Vector of volatilities of two assets
    % rho: Correlation coefficient between the two assets
    % q: Vector of dividend yields of two assets
    % N: Number of time periods
    % H: Barrier level
    % P: Number of paths to simulate

    % Initialize parameters
    mu = r - q - 0.5 * sigma.^2; % Drift term vector
    dt = T / N; % Time step size
    discount_factor = exp(-r * T); % Discount factor for final payoff

    % Generate correlated random numbers for both assets for all paths
    Z1 = randn(P, N);
    Z2 = rho * Z1 + sqrt(1 - rho^2) * randn(P, N);
    
    % Initialize asset price paths
    S1_paths = zeros(P, N);
    S2_paths = zeros(P, N);
    S1_paths(:, 1) = S0(1);
    S2_paths(:, 1) = S0(2);

    % Vectorized path generation
    for t = 2:N
        S1_paths(:, t) = S1_paths(:, t-1) .* exp(mu(1) * dt + sigma(1) * Z1(:, t) * sqrt(dt));
        S2_paths(:, t) = S2_paths(:, t-1) .* exp(mu(2) * dt + sigma(2) * Z2(:, t) * sqrt(dt));
    end

    % Check if the barrier condition was met for each path
    barrier_met = any(S2_paths >= H, 2);

    % Calculate payoffs
    payoff_A = max(S2_paths(:, end) - S1_paths(:, end), 0); % Payoff if barrier is breached
    payoff_A(~barrier_met) = max(X - S1_paths(~barrier_met, end), 0); % Payoff if barrier is not breached

    % Control variate payoff
    CV_final = discount_factor * (0.75 * max(X - S1_paths(:, end), 0) + 0.25 * max(S2_paths(:, end) - X, 0));

    % Monte Carlo estimates
    A = discount_factor * payoff_A;
    estimates = mean(CV_final);

    % Exact control variate using Black-Scholes formula
    exact_CV = 0.75 * BS_put(S0(1), X, r, T, sigma(1), q(1)) + ...
               0.25 * BS_call(S0(2), X, r, T, sigma(2), q(2));

    % Covariance and B calculation
    cov_matrix = cov(A, CV_final);
    B = cov_matrix(1, 2) / var(CV_final);

    % Final variance-reduced estimator
    v = mean(A) - B * (estimates - exact_CV);
end
