function v = FD_eds_call(S0, X, r, T, sigma, q, N, I)
    Smax = 3*X;
    dt = T/N;
    h = Smax/I;
    VGrid = zeros(I+1, N+1);

    % Boundary condition
    VGrid(I+1, :) = (Smax - X)*exp(-r*(T - (0: dt: T)));    % at Smax
    VGrid(1, :) = 0;       %  S = 0

    % Terminal condition
    VGrid(:, N+1) = max((0:I)*h - X, 0);

    i = (1: I-1)';
    i_sq = i.^2;

    % Explicit scheme II
    c = (0.5*sigma^2 * i_sq + 0.5*(r-q)*i) *dt / (1 + r*dt);
    b = (1 - sigma^2 * i_sq * dt) / (1 + r*dt);
    a = (0.5*sigma^2 * i_sq - 0.5*(r-q)*i) *dt / (1 + r*dt);

    % check for negative coeff
    len01 = length(find(a < 0));
    len02 = length(find(b < 0));
    len03 = length(find(c < 0));

    if (len01 > 0 || len02 > 0 || len03 > 0 )
        disp(['# of -ve coeff a: ', num2str(len01)]);
        disp(['# of -ve coeff b: ', num2str(len02)]);
        disp(['# of -ve coeff c: ', num2str(len03)]);
    end

    ishift = 1;

    % backward time 
    for n = N: -1 : 1
        VGrid(i+ishift, n) = a .* VGrid(i-1 + ishift, n+1) ...
            + b .* VGrid(i + ishift, n+1) ...
            + c .* VGrid(i+1 + ishift, n+1);
    end
    
    j = floor(S0/h);
    % disp([num2str(VGrid(:, 1))]);
    v = (S0 - j*h)/h * VGrid(floor(S0/h) + ishift, 1) + ((j+1)*h - S0)/h *VGrid(ceil(S0/h) +ishift, 1);
end