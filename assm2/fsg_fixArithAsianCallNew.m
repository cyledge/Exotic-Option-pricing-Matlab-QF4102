function v = fsg_fixArithAsianCallNew(S0, X, r, T, sigma, q, N, L)
    % eg L = 2, p (rho) = 1/2
    % Precompute constants
    dt = T / N;                 % Time step size
    dx = sigma * sqrt(dt);       % Log step size
    u = exp(sigma * sqrt(dt));   % Up movement factor
    p = (exp((r - q) * dt) - (1/u)) / (u - (1/u));  % Risk-neutral probability

    % Initialize 3D grids (for each time step, stock price, and second state variable)
    % outest loop is t
    % (t, k, i)
    % transforming i to index: ceil((i+n+1)/2)
    % transforming k to index: k + L*n + 1

    ishift =  1;
    kshift = L*N + 1;

    % average at terminal
    A = zeros(2*L*N+1, 1); 
    for k = -L*N: L*N
        A(k + kshift) = S0 * exp(k/L * dx);
    end

   
    % price at terminal
    Vn1 = zeros(2*L*N+1, N+1);  % (k, i)

    for i = 0:N   
        for k = -L*N: L*N
            Vn1(k + kshift, i + ishift) = max(A(k+kshift) - X, 0);
        end
    end


    for n = (N-1):-1:0 % time steps  
        Vn = zeros(2*L*(n-1) + 1, n+1);   % value at this time step

        for i = 0: 1: n     % for different S
            S_n = S0 * u^(2*i - n);

            for k = -L*n: L*n   % for different A
                % Compute average at time step n
                A_n = A(k+ kshift); 
                
                % i+1 is used
                A_n1_ku = (S_n * u + (n+1) * A_n) / (n+2);      % up avg
                A_n1_kd = (S_n * (1/u) + (n+1) * A_n) / (n+2);    % down avg

                k_n1_u = log(A_n1_ku/S0) * L /dx;    % k's are index over [-LN, -LN+1, ..., LN]
                k_n1_d = log(A_n1_kd/S0) * L /dx;
                
                k_fl_u = floor(k_n1_u);   k_cl_u = ceil(k_n1_u);
                k_fl_d = floor(k_n1_d);   k_cl_d = ceil(k_n1_d);
                
                
                % interpolate
                if k_fl_u == k_cl_u
                    ku_index = k_n1_u + kshift; % index of matrix Vn1 (previous V (2*(n+1)*L+1, (n+1)+1))
                    v_n1_ku = Vn1(ku_index, i+1 + ishift);
                else 
                    ku_index = k_fl_u + kshift;                    
                    v_n1_ku_fl = Vn1(ku_index, i+1 +ishift);
                    v_n1_ku_cl = Vn1(ku_index+1, i+1 + ishift);
                    v_n1_ku = (k_n1_u - k_fl_u) * v_n1_ku_fl +  (k_cl_u - k_n1_u) * v_n1_ku_cl;
                end

                if k_fl_d == k_cl_d
                    kd_index = k_n1_d + kshift; % index of matrix Vn1 (previous V (2*(n+1)*L+1, (n+1)+1))
                    v_n1_kd = Vn1(kd_index, i+ ishift);
                else 
                    kd_index = k_fl_d + kshift;
                    v_n1_kd_fl = Vn1(kd_index, i +ishift);
                    v_n1_kd_cl = Vn1(kd_index+1, i +ishift);
                    v_n1_kd = (k_n1_d - k_fl_d) * v_n1_kd_fl +  (k_cl_d - k_n1_d) * v_n1_kd_cl;
                end

                Vn(k + kshift, i + ishift) = exp(-r * dt) * (p * v_n1_ku + (1 - p) * v_n1_kd);
                
            end
        end
        Vn1 = Vn;
    end
    
    % disp(Vn);
    v = Vn(kshift, 1);
end


