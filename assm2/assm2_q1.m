 %% 1(ii)

S0 = 1; X = 1;
T = 0.5;
r = 0.02;
sigma = 0.5;
q = 0.03;
dt = 0.01;
N = T/dt;    % N = 50
h = 0.05;
I = 3*X/h;
v_fd_eds = FD_eds_call(S0, X, r, T, sigma, q, N, I);    % -5.6011e+20
exact_v = BS_call(S0, X, r, T, sigma, q);     % exact 0.1361


%% 1(iii)

for N = 100: 100: 10000
    tem_v = FD_eds_call(S0, X, r, T, sigma, q, N, I);
    if abs(tem_v - exact_v) < 0.1
        break;
    end
    disp(tem_v);
end
disp(N);

%%
v_fd_eds = FD_eds_call(S0, X, r, T, sigma, q, 436, I);   %0.1358
disp(v_fd_eds)
% 436 is the lower bound

%% 1(iv)
% v_fd_eds = 0.1358
% exact_v = 0.1361 > 0.1358 = v_fd_eds
% one of the reasons is due to the Smax, we trancate S to [0, 3*X]

%% 1(v)
v_array = [];
N_array = [];

for N = 436: -5: 0
    tem_v = FD_eds_call(S0, X, r, T, sigma, q, N, I);
    v_array = [v_array, tem_v];
    N_array = [N_array, N];
    if tem_v > 1 || tem_v < 0.05
        disp(tem_v);
        break;
    end
end 
disp(N) % N = 336

for N = 341: -1: 336
    tem_v = FD_eds_call(S0, X, r, T, sigma, q, N, I);
    v_array = [v_array, tem_v];
    N_array = [N_array, N];  
end

figure
hold on
yline(exact_v, 'r--');
plot(N_array, v_array , '.')
legend('BS', 'FD eds')
xlabel('N')
ylabel('price')
title('FD_eds price as N decrease')
hold off

% starting from 338, the value become negative which is meaningless