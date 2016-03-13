%%
clc
clear
close all

%% Polynomial Design
N = 22;

a = [0.0,    0.5];
b = [0.4,    1.0];
m = [0.0,    0.1];

x = pd(N, a, b, m);


%% Plot
n = 0:2*N;
figure, plot(n, real(x), n, imag(x));
xlabel('n');
title('Polynomial');
legend('real','imaginary');

i = 0:2*N;
t = linspace(0,1,1000)';
T = bsxfun(@power, t, i);
X = T * x;

figure, plot(t, X);
xlabel('t');
ylabel('Reponse');
title('Response');
