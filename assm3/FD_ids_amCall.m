function v = FD_ids_amCall(S0, X, r, T, sigma, q, I, N, xmin, xmax, omega, eps)
    dt = T/N;
    dx = (xmax - xmin)/I;
    UGrid = zeros(I+1, N+1);
    sigma_sq = sigma^2;
    % xi_n = xmin + i*dx, i = 0:I

    % boundary condition
    UGrid(I+1, :) = X*(exp(xmax) - 1).*exp(-r.*(T - (T: -dt: 0)));    % at Smax/xmax, len(N+1)
    UGrid(1, :) = X*max(exp(xmin) - 1, 0);          % at xmin, should be = 0 otherwise xmin is meaningless
    
    % terminal condition
    %VGrid(:, N+1) = max((0:I)*h - X, 0); % h = Smax/I;

    UGrid(:, 1) = X * max(exp(xmin + (0:I)*dx) - 1, 0);

    a = 0.5*sigma_sq*dt/(dx^2);
    b = 1 - r*dt - dt/dx * (r - q + 0.5*sigma_sq) - sigma_sq*dt/(dx^2);
    c = dt/dx * (1/(2*dx) * sigma_sq + r - q + 0.5*sigma_sq);

    ishift = 1;
    i = (1: I-1)';
    %i_sq = i.^2;

    % Set up sparse tri-diagonal matrix of coefficients
    CoeffMatrix = spdiags([c, b, a], -1:1, I-1, I-1)';
    
    %xs = xmin + i*dx; % len = I-1
    %payoff_n = X * max(exp(xs) - 1, 0);
    %disp("first grid");
    %disp(UGrid);
    % "backward" which is forward time 
    for n = 2:N+1 % backward time recursive
        % RhsB = V_n+1 + Fn;  len = I-1
        RhsB = UGrid(i+ishift, n-1); % set up right hand side vector
        RhsB(1) = RhsB(1) - a *UGrid(0+ishift, n);
        RhsB(I-1) = RhsB(I-1) - c *UGrid(I+ishift, n);

        xs = xmin + i*dx;   % len = I-1
        payoff_n = X * max(exp(xs) - 1, 0);
        
        Vn1 = CoeffMatrix\RhsB;   % solve linear system
        %payoff_n = exp(-r * (T-n));
        %disp("Vn");
        %disp(Vn);
        UGrid(i+ishift, n) = projSOR(CoeffMatrix, RhsB, Vn1, I - 1, eps, omega, payoff_n);

    end
    % disp("last grid");
    % disp(UGrid);
    
    if log(S0/X) < xmin %if it is really this case, xmin and xmax are not well chosen
        v = UGrid(1, N+1);
    elseif log(S0/X) > xmax
        v = UGrid(I+1, N+1);
    else 
        S0index = (log(S0/X) - xmin)/dx;
        v = (S0index - floor(S0index)) * UGrid(floor(S0index)+ishift, N+1) + (ceil(S0index) - S0index) * UGrid(ceil(S0index)+ishift, N+1);
    end
    
end

%%
function V_n1 = projSOR(A, b, Vn, M, epsilon, omega, payoff_n)
    % pay_off_n + Vn --> Vn1
    % M = len(Vn) = I - 1
    xk = Vn; % Initial guess
    converged = false;
    % k is useless

    while ~converged
        x_k1 = xk; % Prepare next iteration vector
        
        for i = 1:M
            % x^(k+1)_(i,gs)
            sum1 = sum(A(i, 1: i-1) .* x_k1(1: i-1)');    % summation for j = 1: i-1
            sum2 = sum(A(i, i+1: M) .* xk(i+1: M)');       % Summation for j = i+1: M
            x_k_igs = (1/ A(i,i)) * (-sum1 - sum2 + b(i));
            
            % Projected SOR update
            x_k1(i) = max((1 - omega) * xk(i)  +  omega * x_k_igs, payoff_n(i));
        end
        
        % Check for convergence
        if norm(x_k1 - xk, inf) < epsilon
            converged = true;
        end

        xk = x_k1; % Update x_k for the next iteration

    end

    % Set V_n as the converged solution
    V_n1 = xk;
end