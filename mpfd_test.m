%%
clc
clear
close all

%% Minimum Phase Filter Design
N = 22;

a = [-pi/2,     -pi,        2*pi/3];
b = [pi/2,      -2*pi/3,    pi];
l = [0.97,       0.0,        0.0];
u = [1.0,       0.03,       0.03];

[x, G] = mpfd(N, a, b, l, u);

%% Plot
n = 0:N;
figure, plot(n, real(x), n, imag(x));
xlabel('n');
title('Filter');
legend('real','imaginary');

w = linspace(-pi, pi, 1000);
X = fftshift(fft(x, 1000));
figure, plot(w, abs(X));
xlabel('\omega');
ylabel('Magnitude');
title('Magnitude response');
