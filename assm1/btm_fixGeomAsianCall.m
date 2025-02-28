function c = btm_fixGeomAsianCall(S0, X, r, T, sigma, q, N)
    % fix geom asian call
    % max(A - X, 0)
    dt = T/N;
    u = exp(sigma*sqrt(dt));
    d = 1/u;
    p =  ( exp((r-q)*dt)-d ) / (u-d);
    jshift = 1;

    S = zeros(2^N, N + 1);
    S(1, 1) = S0;
    for n = 2:N+1
        for j = 1:(2^(n-2))
            S(2*j - 1, n) = S(j, n-1)*u;
            S(2*j, n) = S(j, n-1)*d;
        end
    end 
    %disp('S');
    %disp(S);
    
    geom = zeros(2^N, N + 1); % store ln(A)
    geom(1, 1) = log(S0);
    for n = 2:N+1
        for j = 1:(2^(n-1))
            %geom(j, n) = ( S(floor(j/2)+1, n) * geom(ceil(j/2), n-1)^(n-1) )^(n);
            geom(j, n) = (log(S(j, n))    +  (n-1)*geom(ceil(j/2), n-1)  )/n;
        end
    end
    %disp(geom);

    % terminal time
    V = max(exp(geom(:, N + 1)) - X, 0);
    %disp(V)

    % backward iterations
    for n = N-1:-1:0
        i = 2^n;
        j = 1:i;
        V(j) = exp(-r*dt)*(p*V(2*j - 1) + (1-p)*V(2*j));
        %disp('for');
        %disp(n);
        %disp(V);
    end
    % output V_0^0
    c = V(1);
    
end
