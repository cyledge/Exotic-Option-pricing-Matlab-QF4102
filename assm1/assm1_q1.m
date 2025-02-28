%% 1(ii)
S0 = 0.8: 0.01: 1.8;
X = 1.2;
sigma = 0.3;
r = 0.05;
T = 1;
q = 0.01;
H = 0.9*ones(1, 101);
c_do = BS_doCall(S0, X, r, T, sigma, q, H);

figure
hold on
plot(S0, c_do, '--')
xlabel('S0')
ylabel('c')
title('down-out barrier call option value, H = 0.9')
hold off

%% 1(iii)
figure
hold on
c = BS_call(S0, X, r, T , sigma, q);
plot(S0, c_do, '--')
plot(S0, c, '.')
legend('down-out', 'vanilla')
xlabel('S0')
ylabel('price')
title('call option values of down-out barrier and vanilla, H = 0.9')
hold off


%% 1(iv)
H = 0.4: 0.01: 1.2;
S0 = 1.3 * ones(1, 81);
c_do = BS_doCall(S0, X, r, T, sigma, q, H);
figure
hold on
plot(H, c_do, '--')
xlabel('H')
ylabel('price')
title('call option values of down-out barrier, S0 = 1.3')
hold off

%%
figure
hold on
c = BS_call(S0, X, r, T , sigma, q);
plot(H, c_do, '--')
plot(H, c, '.')
legend('down-out', 'vanilla')
xlabel('H')
ylabel('price')
title('call option values of down-out barrier and vanilla, S0 = 1.3')
hold off

%% 1(v)
H = 0.9;
S0 = 1.3;
%N = 3070;
c_BTM_do = 3070:3180;
for N = 3070:3180
    c_BTM_do(N - 3069) = btm_doCall(S0, X, r, T, sigma, q, H, N); 
end

c_BS_do = BS_doCall(S0, X, r, T, sigma, q, H);

N = 3070:3180;
figure
hold on
plot(N, c_BTM_do - c_BS_do, '--')
xlabel('N')
ylabel('error')
title('BTM model error compared to BS price')
hold off

%% 1(vi)

% smallest error 
[error, N] = min(abs(c_BTM_do - c_BS_do));
N = 3070 + N - 1;  % N = 3077
disp('first local minima')
disp(N);

% another local minima, 2nd smallest error
[error, N] = min(abs(c_BTM_do(80:111) - c_BS_do));
N = 3070 + 80 + N - 1;  % N = 3166
disp('second local minima');
disp(N);

% first jump
[error, N] = min(c_BTM_do - c_BS_do);
N = 3070 + N - 1;  % N = 3167
disp('first jump');
disp(N);

% another local minima, 2nd smallest error
[error, N] = min(c_BTM_do(1:80) - c_BS_do);
N = 3070 + N - 1;  % N = 3077
disp('second jump');
disp(N)

%{
N < T(i*sigma)^2 / (ln(H/S0))^2
3167 < (0.3i)^2 / 0.135221513
428.2465317 < (0.3i)^2
--> i < -68.98039 or i > 68.98 (rej.)
--> i = -68

3077 < (0.3i)^2 / 0.135221513
--> i < -67.99318557
--> i = -67

%}
