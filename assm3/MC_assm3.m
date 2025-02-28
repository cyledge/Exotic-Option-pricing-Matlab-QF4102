function v = MC_assm3(S0, X, r, T, sigma, rho, q, N, H, P)
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

    % Matrix to store final asset prices for all paths
    price_sum = 0;

    % Generate correlated random variables for both assets
    for i = 1:P
        marker = 0; % Reset the marker for each path
        
        % Initialize asset prices for the current path
        S1 = S0(1);
        S2 = S0(2);

        % Simulate the path for both assets
        for t = 1:N
            % Generate correlated random numbers
            random1 = randn;
            random2 = rho * random1 + sqrt(1 - rho^2) * randn;
            
            % Update asset prices using geometric Brownian motion
            S1 = S1 * exp(mu(1) * dt + sigma(1) * random1 * sqrt(dt));
            S2 = S2 * exp(mu(2) * dt + sigma(2) * random2 * sqrt(dt));
            
            % Check if the barrier condition is met
            if S2 >= H
                marker = 1;
            end
        end

        % Calculate the payoff based on whether the barrier was breached
        if marker == 1
            payoff = max(S2 - S1, 0);
        else
            payoff = max(X - S1, 0);
        end

        % Accumulate the discounted payoff
        price_sum = price_sum + payoff;
    end

    % Calculate the average payoff and discount it
    v = discount_factor * (price_sum / P);
end
