% Parameters
S0 = 0.85;       % Current underlying price
r = 0.03;        % Risk-free rate
T = 0.5;         % Time to maturity
sigma = 0.35;    % Volatility
q = 0;           % Dividend yield
% runmin = 0.75;   % Running minimum 

N_values = 200:200:10000;  % Range of N from 200 to 10000 in increments of 200
option_values = zeros(size(N_values));  

% Loop over different values of N
for idx = 1:length(N_values)
    N = N_values(idx);
    option_values(idx) = btm_1vFlLookbackCall(S0, r, T, sigma, q, N, runmin);
end

% Plotting the option values versus N
figure;
plot(N_values, option_values, 'LineWidth', 2);
xlabel('Number of Time Steps (N)');
ylabel('Option Value');
title('Option Value vs. Number of Time Steps (N) for Lookback Call');
grid on;

function v = btm_1vFlLookbackCall(S0, r, T, sigma, q, N, runmin) 

% Check if runmin is not provided, set it to S0
    if nargin < 7
        runmin = S0;  % Default value if runmin is not specified
    end

    dt = T / N;                        
    dx = sigma * sqrt(dt);    % Step size 
    u = exp(dx);                       
    d = 1 / u;                        
    df = exp(-r * dt);               
    p = (exp((r - q) * dt) - d) / (u - d);   
    x0 = log(min(S0, runmin) / S0);  
    j = max(0, floor(x0 / dx));       
 
    % Initialization 
    jshift = 1; 
    nshift = 1; 
    initial = max(j - N, 0); 
    V = zeros(j + N + 2, N + 1);  
     
    % Terminal condition 
    i = initial:(j + N + 1); 
    S = S0 .* exp(i .* dx);         
    V(i + jshift, N + 1) = max(S - runmin, 0); 
     
    % Backward recursive process 
    for n = N-1:-1:0 
        if (j - n <= 0) 
            V(1, n + nshift) = df * (p*u * V(1 + jshift, n+1 + nshift) + (1-p)*d * V(1, n+1 + nshift)); 
        end 
        i = max((j - n), 1):(j + n + 1);  
        V(i+jshift, n + nshift) = df * (p*u * V(i+1 + jshift, n+1 + nshift) + (1-p)*d * V(i-1 + jshift, n+1 + nshift)); 
    end 
 
    % Interpolation for the final option value 
    y1 = S0 * V(j + jshift, nshift); 
    y2 = S0 * V(j + 1 + jshift, nshift); 
    x1 = j * dx; 
    x2 = (j + 1) * dx; 
    v = ((x2 - x0) / dx) * y1 + ((x0 - x1) / dx) * y2; 
end
