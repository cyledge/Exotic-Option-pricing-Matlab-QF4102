function c = btm_doCall(S0, X, r, T, sigma, q, H, N) 
    dt = T/N;
    u = exp(sigma*sqrt(dt));
    % disp('u');
    % disp(u);
    d = 1/u;
    p =  (exp((r-q)*dt)-d ) / (u-d);

    jshift = 1;

    % terminal time
    j = 0:N;
    S = S0 * u.^(2*j-N); % dim = (1, N+1)
    % disp('S');
    % disp(S);
    % H <= X --> S <= H <= X, 
    V = max(S-X, 0); % dim = (1, N+1)

    % backward iterations
    for n = N-1:-1:0
        j = 0:n;
        % V change from (1, n+2) --> (1, n+1)
        V(j+jshift) = exp(-r*dt)*( p*V(j+jshift+1) + (1-p)*V(j+jshift) );
        % disp(V);
        
        % every option with index < smallj will be voided, as St < H
        smallj = (log(H/S0)/log(u) + n)/2 ;
        if smallj > 0
            smallj = floor(smallj) - (floor(smallj) == smallj) + 1;
            V(1:smallj) = 0;
           
        end
       
        %{
        S_tem = S0 * u .^(2*j - n);
        V(S_tem < H) = 0;
        disp('S_tem');
        disp(S_tem);
        %}

    end
    % output V_0^0
    c = V(1);
end