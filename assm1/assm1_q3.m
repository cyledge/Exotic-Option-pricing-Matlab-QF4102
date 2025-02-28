S0 = 95;
X = 90;
r = 0.04;
T = 0.5;
sigma = 0.3;
q = 0;

c_fixAsian = zeros(1, 4);
time = zeros(1, 4);
i = 1;
for N = 5:5:20
    tic;
    c_fixAsian(i) = btm_fixGeomAsianCall(S0, X, r, T, sigma, q, N);
    time(i) = toc;
    i = i + 1;
    
end

N = 5:5:20;
figure
hold on
plot(N, c_fixAsian, '--')
xlabel('N')
ylabel('price')
title('BTM fixGeomAsianCall vs N')
hold off

figure
hold on
plot(N, time, '--')
xlabel('N')
ylabel('time')
title('time vs N')
hold off