function c = btm_EurCall(S0,X,r,T,sigma,q,N)
    dt = T/N;
    u = exp(sigma*sqrt(dt));
    d = 1/u;
    p =  ( exp((r-q)*dt)-d ) / (u-d);

    jshift = 1;

    % terminal time
    j = 0:N;
    S = S0 * u.^(2*j-N);
    V = max(S-X, 0);

    % backward iterations
    for n = N-1:-1:0
        j = 0:n;
        V(j+jshift) = exp(-r*dt)*( p*V(j+jshift+1) + (1-p)*V(j+jshift) );
        disp(n);
        disp(V);
    end
    % output V_0^0
    c = V(1);
end